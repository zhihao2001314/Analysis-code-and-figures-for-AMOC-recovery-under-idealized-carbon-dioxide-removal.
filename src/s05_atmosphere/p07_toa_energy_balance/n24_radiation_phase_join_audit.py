# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import netCDF4 as nc
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
stage = root / 'output/s05_atmosphere/p07_toa_energy_balance'
destination = stage / 'n24_radiation_phase_join_audit'
assert not destination.exists(), '正式结果存在，须先归档。'
models = ['CESM2', 'CNRM-ESM2-1', 'CAS-ESM2-0', 'GFDL-ESM4', 'NorESM2-LM']
quantities = ['rsdt', 'rsut', 'rlut', 'asr', 'negative_olr', 'eei']
starts = {'CESM2': [1, 1], 'CNRM-ESM2-1': [1850, 1990], 'CAS-ESM2-0': [200, 340],
          'GFDL-ESM4': [1, 141], 'NorESM2-LM': [1, 141]}
old = pd.concat([pd.read_csv(stage / 'n11_picontrol_corrected_eei/corrected_annual.csv', float_precision='round_trip'),
                 pd.read_csv(stage / 'n21_cas_picontrol_corrected_radiation/corrected_annual.csv', float_precision='round_trip')], ignore_index=True)
gfdl = pd.read_csv(stage / 'n11_picontrol_corrected_eei/gfdl_source1_140_sensitivity.csv', float_precision='round_trip')
selected_columns = ['model', 'phase', 'source_year', 'simulation_year', 'control_year']
for quantity in quantities:
    selected_columns += [quantity, quantity + '_control_fit', quantity + '_relative_control']
corrected = pd.concat([old.loc[~old.model.eq('GFDL-ESM4'), selected_columns], gfdl[selected_columns]], ignore_index=True)
monthly = pd.concat([pd.read_csv(stage / 'n01_toa_energy_balance/monthly_global_components_long.csv'),
                     pd.read_csv(stage / 'n20_cas_toa_radiation/monthly_global_long.csv')], ignore_index=True)
inventory = pd.read_csv(stage / 'n01_toa_energy_balance/source_inventory.csv')
cas_inventory = pd.read_csv(stage / 'n20_cas_toa_radiation/source_inventory.csv')
audit_rows, monthly_checks, branch_headers = [], [], []


# %% 3. 每变量140年逐月完整性，独立复核年平均与真实源分支

for model in models:
    model_series = corrected.loc[corrected.model.eq(model)].sort_values('simulation_year')
    np.testing.assert_array_equal(model_series.simulation_year, np.arange(1, 281))
    np.testing.assert_array_equal(np.diff(model_series.control_year), np.ones(279))
    for phase_index, phase in enumerate(['increase', 'removal']):
        first_year = starts[model][phase_index]
        for quantity in ['rsdt', 'rsut', 'rlut']:
            selected = monthly.loc[monthly.model.eq(model) & monthly.phase.eq(phase)
                                    & monthly.variable.eq(quantity) & monthly.source_year.between(first_year, first_year + 139)].copy()
            selected = selected.sort_values(['source_year', 'month'])
            np.testing.assert_array_equal(selected.source_year.to_numpy() * 12 + selected.month.to_numpy() - 1,
                                           np.arange(first_year * 12, (first_year + 140) * 12))
            assert len(selected) == 1680
            selected['weighted'] = selected.global_Wm2 * selected.days
            yearly = selected.groupby('source_year').weighted.sum() / selected.groupby('source_year').days.sum()
            reference = model_series.loc[model_series.phase.eq(phase)].set_index('source_year')[quantity]
            error = float(np.max(np.abs(yearly.to_numpy() - reference.to_numpy())))
            assert error < 1e-10
            monthly_checks.append({'model': model, 'phase': phase, 'variable': quantity,
                                   'source_first': first_year, 'source_last': first_year + 139,
                                   'months': len(selected), 'max_annual_error_Wm2': error})
        if model == 'CAS-ESM2-0':
            path = Path(cas_inventory.loc[cas_inventory.experiment.eq(['1pctCO2', '1pctCO2-cdr'][phase_index])
                                          & cas_inventory.variable.eq('rsdt'), 'path'].iloc[0])
        else:
            path = Path(inventory.loc[inventory.model.eq(model) & inventory.variable.eq('rsdt')
                                      & inventory.experiment.eq(['1pctCO2', '1pctCO2-cdr'][phase_index]), 'path'].iloc[0])
        with nc.Dataset(path) as dataset:
            time = dataset['time']
            parent_date = nc.num2date(dataset.branch_time_in_parent, dataset.parent_time_units, calendar=time.calendar)
            child_date = nc.num2date(dataset.branch_time_in_child, time.units, calendar=time.calendar)
            branch_headers.append({'model': model, 'phase': phase, 'path': str(path),
                                   'dimensions': list(dataset['rsdt'].dimensions),
                                   'time_units': time.units, 'calendar': time.calendar,
                                   'parent_experiment': dataset.parent_experiment_id,
                                   'parent_units': dataset.parent_time_units,
                                   'parent_branch_date_literal': str(parent_date),
                                   'child_branch_date_literal': str(child_date),
                                   'branch_method': dataset.branch_method})
            if model == 'GFDL-ESM4' and phase == 'removal':
                assert parent_date.year == child_date.year == 141
    print(model, '280 years monthly coverage and annual means verified', flush=True)


# %% 4. 同一280年序列连续11年平滑，不移动原始值或截距

for model in models:
    selected = corrected.loc[corrected.model.eq(model)].sort_values('simulation_year').copy()
    previous = old.loc[old.model.eq(model)].sort_values('simulation_year').set_index('simulation_year')
    for quantity in quantities:
        column = quantity + '_relative_control'
        selected[column + '_smooth11'] = selected[column].rolling(11, center=True, min_periods=6).mean()
        new = selected.set_index('simulation_year')
        first, second = new.loc[140], new.loc[141]
        audit_rows.append({'model': model, 'quantity': quantity,
                           'old_source140': int(previous.loc[140, 'source_year']),
                           'new_source140': int(first.source_year), 'new_source141': int(second.source_year),
                           'raw140_Wm2': first[quantity], 'raw141_Wm2': second[quantity],
                           'raw_step_Wm2': second[quantity] - first[quantity],
                           'control_fit_step_Wm2': second[quantity + '_control_fit'] - first[quantity + '_control_fit'],
                           'corrected_step_Wm2': second[column] - first[column],
                           'old_phase_smooth_step_Wm2': previous.loc[141, column + '_smooth11'] - previous.loc[140, column + '_smooth11'],
                           'new_continuous_smooth_step_Wm2': second[column + '_smooth11'] - first[column + '_smooth11']})
        corrected.loc[selected.index, column + '_smooth11'] = selected[column + '_smooth11']
corrected = corrected.sort_values(['model', 'simulation_year']).reset_index(drop=True)
for suffix in ['', '_control_fit', '_relative_control', '_relative_control_smooth11']:
    np.testing.assert_allclose(corrected['eei' + suffix],
                               corrected['asr' + suffix] + corrected['negative_olr' + suffix], atol=1e-10, rtol=0)
unchanged = corrected.loc[~corrected.model.eq('GFDL-ESM4')].sort_values(['model', 'simulation_year'])
original = old.loc[~old.model.eq('GFDL-ESM4')].sort_values(['model', 'simulation_year'])
pd.testing.assert_frame_equal(unchanged[selected_columns].reset_index(drop=True),
                              original[selected_columns].reset_index(drop=True), check_exact=True)
for quantity in quantities:
    np.testing.assert_array_equal(corrected.loc[corrected.model.eq('GFDL-ESM4'), quantity + '_control_fit'].to_numpy(),
                                   gfdl[quantity + '_control_fit'].to_numpy())


# %% 5. 保存断接原因与复核证据

destination.mkdir(parents=True)
corrected.to_csv(destination / 'corrected_annual.csv', index=False)
pd.DataFrame(audit_rows).to_csv(destination / 'seam_audit.csv', index=False)
pd.DataFrame(monthly_checks).to_csv(destination / 'monthly_join_qc.csv', index=False)
(destination / 'branch_headers.json').write_text(json.dumps(branch_headers, ensure_ascii=False, indent=2), encoding='utf-8')
metadata = {'status': 'passed_processing_fixes_CAS_branch_unverified', 'models': models,
            'starts_increase_removal': starts, 'raw_dimensions': 'All rsdt/rsut/rlut(time,lat,lon), native 1Dlat/lon; units W m-2.',
            'causes': ['Plot previously used two separate line artists for each component and phase, leaving an artificial gap.',
                       'Independent phase smoothing switched from years135-140 to141-146 at the seam, creating endpoint steps.',
                       'GFDL radiation inherited source6-145 rise followed by141-280 decline, despite raw radiation availability1-140 and branch date141.'],
            'fix': 'GFDL uses validated source1-140 rise and141-280 decline, unchanged existing control OLS fit. All five radiative series use one 11-year centered window spanning280years, min6 at outer endpoints only; plot one line per280year series.',
            'piControl': 'Same fitted intercepts/slopes; continuous matched control years for every series. CAS mapping remains provisional as before.',
            'no_offset_or_interpolation': True, 'no_missing_months': True,
            'GFDL_scope': 'Only Figure3 radiation corrected here; previous ocean-start6 products are not silently changed.',
            'CAS_limitation': 'Raw phase-year339 to340 global LW contribution changes -0.767Wm-2. No missing year or variable mismatch; <=0.0007Wm-2 control-fit seam change. Parent/child branch-date attributes conflict with actual coordinates, so physical restart continuity remains unverified. Raw fluctuation retained, not shifted away.',
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in
                             [stage / 'n11_picontrol_corrected_eei/corrected_annual.csv',
                              stage / 'n11_picontrol_corrected_eei/gfdl_source1_140_sensitivity.csv',
                              stage / 'n21_cas_picontrol_corrected_radiation/corrected_annual.csv',
                              stage / 'n01_toa_energy_balance/monthly_global_components_long.csv',
                              stage / 'n20_cas_toa_radiation/monthly_global_long.csv']]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(audit_rows).loc[pd.DataFrame(audit_rows).quantity.isin(['eei', 'asr', 'negative_olr'])].to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
