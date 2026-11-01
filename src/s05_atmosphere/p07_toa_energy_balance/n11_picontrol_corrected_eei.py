# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import netCDF4 as nc
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
stage = Path('s05_atmosphere/p07_toa_energy_balance')
control_directory = root / 'output' / stage / 'n10_picontrol_toa'
forced_directory = root / 'output' / stage / 'n01_toa_energy_balance'
destination = root / 'output' / stage / 'n11_picontrol_corrected_eei'
assert not destination.exists(), '正式结果存在，须先归档。'
models = ['CESM2', 'CNRM-ESM2-1', 'NorESM2-LM', 'GFDL-ESM4']
quantities = ['rsdt', 'rsut', 'rlut', 'asr', 'negative_olr', 'eei']
windows = {'initial': [1, 10], 'peak': [131, 140], 'early_removal': [141, 150], 'late': [271, 280]}
forced_metadata = json.loads((forced_directory / 'metadata.json').read_text(encoding='utf-8'))
forced_all = pd.read_csv(forced_directory / 'annual_global_energy.csv', float_precision='round_trip')
forced_inventory = pd.read_csv(forced_directory / 'source_inventory.csv')
gfdl_extra = pd.read_csv(forced_directory / 'gfdl_source_year_sensitivity.csv', float_precision='round_trip')
control_metadata = {}
for model in models:
    control_metadata[model] = json.loads((control_directory / model / 'metadata.json').read_text(encoding='utf-8'))
    assert control_metadata[model]['status'] == 'passed'
fit_rows, branch_rows, control_parts, corrected_parts = [], [], [], []
endpoint_rows, sensitivity_rows, gfdl_parts = [], [], []


# %% 3. 核对父链并拟合对应280年控制背景

for model in models:
    control = pd.read_csv(control_directory / model / 'annual_control.csv', float_precision='round_trip')
    fit_first, fit_last = control_metadata[model]['fit_window']
    selected = control[control.source_year.between(fit_first, fit_last)].copy()
    assert len(selected) == 280 and selected[quantities].notna().all().all()
    np.testing.assert_array_equal(selected.source_year, np.arange(fit_first, fit_last + 1))
    origin = (fit_first + fit_last) / 2
    design = np.column_stack([np.ones(len(selected)), selected.source_year.to_numpy() - origin])
    coefficients = {}
    for quantity in quantities:
        coefficient = np.linalg.lstsq(design, selected[quantity].to_numpy(), rcond=None)[0]
        coefficients[quantity] = coefficient
        control[quantity + '_fit'] = coefficient[0] + coefficient[1] * (control.source_year - origin)
        fit_rows.append({'model': model, 'quantity': quantity, 'first_year': fit_first,
                         'last_year': fit_last, 'n_years': len(selected), 'origin_year': origin,
                         'intercept_at_origin_Wm2': coefficient[0], 'slope_Wm2_per_year': coefficient[1],
                         'slope_Wm2_per_century': coefficient[1] * 100,
                         'control_mean_Wm2': selected[quantity].mean(),
                         'residual_sample_SD_Wm2': np.std(selected[quantity] - design @ coefficient, ddof=1)})
    np.testing.assert_allclose(coefficients['eei'], coefficients['rsdt'] - coefficients['rsut'] - coefficients['rlut'], atol=1e-11, rtol=0)
    control_parts.append(control)

    phase_alignment = {}
    for phase, experiment in [('increase', '1pctCO2'), ('removal', '1pctCO2-cdr')]:
        source_row = forced_inventory[(forced_inventory.model == model) &
                                      (forced_inventory.experiment == experiment) &
                                      (forced_inventory.variable == 'rsdt')].iloc[0]
        with nc.Dataset(source_row.path) as dataset:
            time = dataset.variables['time']
            assert time.calendar == control_metadata[model]['source_files'][0]['calendar']
            parent_branch = nc.num2date(dataset.branch_time_in_parent, dataset.parent_time_units,
                                       calendar=time.calendar)
            child_declared = nc.num2date(dataset.branch_time_in_child, time.units, calendar=time.calendar)
            actual_first = nc.num2date(dataset.variables[time.bounds][0, 0], time.units, calendar=time.calendar)
            assert parent_branch.month == child_declared.month == actual_first.month == 1
            assert dataset.parent_variant_label == control_metadata[model]['member']
            expected_parent = 'piControl' if phase == 'increase' else '1pctCO2'
            assert dataset.parent_experiment_id == expected_parent
            if phase == 'increase':
                assert parent_branch.year == fit_first
                increase_child_origin = child_declared.year
                increase_control_origin = parent_branch.year
                phase_alignment[phase] = [increase_control_origin, increase_child_origin]
                child_used = child_declared.year
                alignment_note = 'Decoded parent and child branch dates; no override.'
            else:
                child_used = child_declared.year
                alignment_note = 'Decoded parent and child branch dates; no override.'
                if model == 'NorESM2-LM':
                    # 原属性child=0解码为0001，但branch_method、父日期及真实首界均为0141。
                    assert child_declared.year == 1 and actual_first.year == parent_branch.year == 141
                    assert '0141-01-01' in dataset.branch_method
                    child_used = actual_first.year
                    alignment_note = 'Explicit interpretation of inconsistent child attribute: use actual0141first bound, consistent with parent date and branch_method. Raw metadata unchanged.'
                assert child_used == forced_metadata['starts_increase_removal'][model][1]
                removal_control_origin = increase_control_origin + parent_branch.year - increase_child_origin
                phase_alignment[phase] = [removal_control_origin, child_used]
            branch_rows.append({'model': model, 'phase': phase, 'path': source_row.path,
                                'parent_experiment_id': dataset.parent_experiment_id,
                                'parent_time_units': dataset.parent_time_units,
                                'branch_time_in_parent': float(dataset.branch_time_in_parent),
                                'branch_time_in_child': float(dataset.branch_time_in_child),
                                'time_units': time.units, 'calendar': time.calendar,
                                'parent_branch_date': str(parent_branch), 'child_declared_date': str(child_declared),
                                'actual_first_date': str(actual_first), 'child_origin_used': child_used,
                                'control_origin_used': phase_alignment[phase][0],
                                'branch_method': dataset.branch_method, 'alignment_note': alignment_note})

    # 预先定义窗口敏感性，不根据结果挑选最有利控制段。
    common = control.dropna(subset=quantities)
    fit_options = {'matched280': selected, 'matched_first140': selected.iloc[:140],
                   'matched_last140': selected.iloc[-140:], 'all_common_available': common}
    main_forced = forced_all[forced_all.model == model].copy()
    forced_control_years = np.empty(len(main_forced), dtype=int)
    for phase in ['increase', 'removal']:
        mask = (main_forced.phase == phase).to_numpy()
        control_origin, child_origin = phase_alignment[phase]
        forced_control_years[mask] = control_origin + main_forced.loc[mask, 'source_year'].to_numpy() - child_origin
    early_control_year = forced_control_years[main_forced.simulation_year.between(1, 10)].mean()
    late_control_year = forced_control_years[main_forced.simulation_year.between(271, 280)].mean()
    raw_delta = main_forced.loc[main_forced.simulation_year.between(271, 280), 'eei'].mean() - main_forced.loc[main_forced.simulation_year.between(1, 10), 'eei'].mean()
    for option, fit_part in fit_options.items():
        x = fit_part.source_year.to_numpy() - origin
        option_coefficient = np.linalg.lstsq(np.column_stack([np.ones(len(x)), x]), fit_part.eei.to_numpy(), rcond=None)[0]
        drift_delta = option_coefficient[1] * (late_control_year - early_control_year)
        sensitivity_rows.append({'model': model, 'fit_option': option, 'fit_first_year': int(fit_part.source_year.min()),
                                 'fit_last_year': int(fit_part.source_year.max()), 'n_years': len(fit_part),
                                 'slope_Wm2_per_century': option_coefficient[1] * 100,
                                 'forced_raw_late_minus_early': raw_delta,
                                 'control_fitted_late_minus_early': drift_delta,
                                 'corrected_late_minus_early': raw_delta - drift_delta,
                                 'control_year_separation': late_control_year - early_control_year})


    # 分别保存仅去漂移与相对于控制拟合的序列。

    cases = ['existing_window', 'gfdl_source1_140'] if model == 'GFDL-ESM4' else ['existing_window']
    for case in cases:
        keep = ['model', 'phase', 'source_year', 'simulation_year'] + quantities
        forced = main_forced[keep].copy()
        if case == 'gfdl_source1_140':
            initial_phase = gfdl_extra[gfdl_extra.source_year.between(1, 140)][keep].copy()
            initial_phase['simulation_year'] = initial_phase.source_year
            forced = pd.concat([initial_phase, forced[forced.phase == 'removal']], ignore_index=True)
        forced = forced.sort_values('simulation_year').reset_index(drop=True)
        assert len(forced) == 280 and forced[quantities].notna().all().all()
        forced['case'] = case
        forced['control_year'] = 0
        for phase in ['increase', 'removal']:
            mask = forced.phase == phase
            control_origin, child_origin = phase_alignment[phase]
            forced.loc[mask, 'control_year'] = control_origin + forced.loc[mask, 'source_year'] - child_origin
        assert forced.control_year.between(fit_first, fit_last).all()
        reference_control_year = forced.loc[forced.simulation_year.between(1, 10), 'control_year'].mean()
        # 阶段内部完整11年，端点只用该阶段可用的6–10年；不跨转折平滑。
        for phase in ['increase', 'removal']:
            mask = forced.phase == phase
            forced.loc[mask, 'smoothing_window_years'] = forced.loc[mask, 'eei'].rolling(11, center=True, min_periods=6).count()
        forced['smoothing_window_years'] = forced.smoothing_window_years.astype(int)
        for quantity in quantities:
            intercept, slope = coefficients[quantity]
            forced[quantity + '_control_fit'] = intercept + slope * (forced.control_year - origin)
            # 前者去掉控制拟合的截距及趋势；后者只减相对初始窗口的趋势，保留绝对偏移。
            forced[quantity + '_relative_control'] = forced[quantity] - forced[quantity + '_control_fit']
            forced[quantity + '_trend_only'] = forced[quantity] - slope * (forced.control_year - reference_control_year)
            for suffix in ['', '_relative_control', '_trend_only']:
                name = quantity + suffix
                forced[name + '_initial_anomaly'] = forced[name] - forced.loc[forced.simulation_year.between(1, 10), name].mean()
                for phase in ['increase', 'removal']:
                    mask = forced.phase == phase
                    forced.loc[mask, name + '_smooth11'] = forced.loc[mask, name].rolling(11, center=True, min_periods=6).mean()
        np.testing.assert_allclose(forced.eei_relative_control,
                                   forced.rsdt_relative_control - forced.rsut_relative_control - forced.rlut_relative_control,
                                   atol=1e-10, rtol=0)
        np.testing.assert_allclose(forced.eei_relative_control_initial_anomaly,
                                   forced.eei_trend_only_initial_anomaly, atol=1e-10, rtol=0)
        for label, (first, last) in windows.items():
            part = forced[forced.simulation_year.between(first, last)]
            record = {'model': model, 'case': case, 'window': label, 'first_simulation_year': first,
                      'last_simulation_year': last, 'control_first_year': int(part.control_year.min()),
                      'control_last_year': int(part.control_year.max())}
            for quantity in quantities:
                for suffix in ['', '_control_fit', '_relative_control', '_trend_only']:
                    record[quantity + suffix] = part[quantity + suffix].mean()
            endpoint_rows.append(record)
        if case == 'existing_window':
            corrected_parts.append(forced)
        else:
            gfdl_parts.append(forced)


# %% 4. 保存结果、差值和完整方法记录

destination.mkdir(parents=True)
corrected = pd.concat(corrected_parts, ignore_index=True)
endpoints = pd.DataFrame(endpoint_rows)
contrasts = []
for model in models:
    for case in endpoints.loc[endpoints.model == model, 'case'].unique():
        part = endpoints[(endpoints.model == model) & (endpoints.case == case)].set_index('window')
        record = {'model': model, 'case': case}
        for quantity in quantities:
            for suffix in ['', '_control_fit', '_relative_control', '_trend_only']:
                name = quantity + suffix
                record[name + '_late_minus_early'] = part.loc['late', name] - part.loc['initial', name]
        contrasts.append(record)
pd.DataFrame(fit_rows).to_csv(destination / 'control_fit_coefficients.csv', index=False)
pd.DataFrame(branch_rows).to_csv(destination / 'branch_alignment.csv', index=False)
pd.concat(control_parts, ignore_index=True).to_csv(destination / 'control_annual_with_fit.csv', index=False)
corrected.to_csv(destination / 'corrected_annual.csv', index=False)
endpoints.to_csv(destination / 'endpoint_comparison.csv', index=False)
pd.DataFrame(contrasts).to_csv(destination / 'late_minus_early.csv', index=False)
pd.DataFrame(sensitivity_rows).to_csv(destination / 'fit_window_sensitivity.csv', index=False)
pd.concat(gfdl_parts, ignore_index=True).to_csv(destination / 'gfdl_source1_140_sensitivity.csv', index=False)
source_paths = [forced_directory / 'annual_global_energy.csv', forced_directory / 'gfdl_source_year_sensitivity.csv']
for model in models:
    source_paths.extend([control_directory / model / 'annual_control.csv', control_directory / model / 'metadata.json'])
metadata = {'status': 'passed', 'models': models, 'fit': 'OLS on raw annual values in branch-aligned280year control window; equal annual weights; intercept+slope; no forced-series detrending.',
            'relative_control_definition': 'N_forced(t)-Nhat_control(control_year(t)); removes control offset and linear drift. This is a control-relative flux anomaly, not absolute Earth energy uptake.',
            'trend_only_definition': 'N_forced(t)-b*(control_year(t)-mean_control_year_of_forced_initial10); retains initial absolute baseline.',
            'endpoint_definition': 'late minus early, corrected response=forced difference minus fitted control difference; identical for both correction conventions.',
            'fit_windows': {model: control_metadata[model]['fit_window'] for model in models},
            'smoothing': 'Centered11year within each phase, shortened to6-10available years at phase endpoints; no cross-phase smoothing or imputation. Actual window lengths saved. Display only; endpoint statistics raw annual.',
            'Nor_alignment': 'Child branch_time=0 inconsistent with actual0141coordinate. Parent date, branch_method and actual first bound all indicate0141. Explicit interpretation recorded; original attributes retained unchanged.',
            'GFDL_alignment': 'Primary preserves source6-145 increase and141-280 removal. Control source106-245 then241-380; phase overlap explicit, no claim of uninterrupted280physical years. Separate source1-140 sensitivity saved.',
            'fit_sensitivity': 'Matched280, first140, last140, all available common control years; diagnostic, not statistical confidence intervals; no significance claim.',
            'original_results_unchanged': True,
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths],
            'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(fit_rows).query('quantity == "eei"').to_string(index=False), flush=True)
print(pd.DataFrame(contrasts)[['model', 'case', 'eei_late_minus_early', 'eei_control_fit_late_minus_early',
                              'eei_relative_control_late_minus_early']].to_string(index=False), flush=True)
