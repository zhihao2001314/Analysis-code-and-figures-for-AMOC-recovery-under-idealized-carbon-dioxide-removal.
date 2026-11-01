# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s05_atmosphere/p07_toa_energy_balance/n18_five_model_figure3'
assert not destination.exists() or not any(destination.iterdir()), '正式结果存在，须先归档。'
radiation_source = root / 'output/s05_atmosphere/p07_toa_energy_balance/n11_picontrol_corrected_eei'
cas_source = root / 'output/s05_atmosphere/p07_toa_energy_balance/n20_cas_toa_radiation'
cas_corrected_source = root / 'output/s05_atmosphere/p07_toa_energy_balance/n21_cas_picontrol_corrected_radiation'
join_source = root / 'output/s05_atmosphere/p07_toa_energy_balance/n24_radiation_phase_join_audit'
temperature_source = root / 'output/s03_tem/n13_six_model_global_temperature_groups'
groups = {'non_recovery': ['CESM2', 'CNRM-ESM2-1'],
          'recovery': ['CAS-ESM2-0', 'GFDL-ESM4', 'NorESM2-LM']}
temperature_groups = {'non_recovery': ['CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL'],
                      'recovery': ['CAS-ESM2-0', 'GFDL-ESM4', 'NorESM2-LM']}
models = groups['non_recovery'] + groups['recovery']
fields = ['eei_relative_control', 'asr_relative_control', 'negative_olr_relative_control']
source_files = [radiation_source / 'corrected_annual.csv', radiation_source / 'metadata.json',
                join_source / 'corrected_annual.csv', join_source / 'metadata.json', join_source / 'seam_audit.csv',
                cas_source / 'annual_global_energy.csv', cas_source / 'metadata.json',
                cas_corrected_source / 'corrected_annual.csv', cas_corrected_source / 'metadata.json',
                cas_corrected_source / 'alignment_sensitivity_maximum.csv',
                temperature_source / 'global_temperature_groups.npz', temperature_source / 'metadata.json',
                temperature_source / 'coastline_global_50m.geojson']
membership_rows, group_rows, endpoint_rows = [], [], []


# %% 3. 读取已核对的温度场，按3对3平均，温度不恢复组保留UKESM

temperature_metadata = json.loads((temperature_source / 'metadata.json').read_text(encoding='utf-8'))
with np.load(temperature_source / 'global_temperature_groups.npz') as saved:
    latitude = saved['latitude'].copy()
    longitude = saved['longitude'].copy()
    source_temperature = saved['model_fields'].copy()
temperature_fields = []
for group, members in temperature_groups.items():
    indices = [temperature_metadata['models'].index(model) for model in members]
    selected = source_temperature[indices]
    assert np.isfinite(selected).all()
    temperature_fields.append(selected.mean(axis=0))
    for model in members:
        membership_rows.append({'quantity': 'temperature', 'group': group,
                                'model': model, 'weight': 1 / len(members)})
for group, members in groups.items():
    for model in members:
        membership_rows.append({'quantity': 'radiation', 'group': group,
                                'model': model, 'weight': 1 / len(members)})
temperature_panels = np.stack(temperature_fields + [temperature_fields[0] - temperature_fields[1]])
latitude_bounds = np.deg2rad(np.column_stack([latitude - .5, latitude + .5]))
weights = np.broadcast_to((np.sin(latitude_bounds[:, 1]) - np.sin(latitude_bounds[:, 0]))[:, None],
                          (180, 360)).copy()
weights /= weights.sum()
for name, values in zip(['non_recovery', 'recovery', 'non_recovery_minus_recovery'], temperature_panels):
    endpoint_rows.append({'quantity': 'tas_residual_degC', 'series': name,
                          'value': float(np.sum(values * weights))})
temperature_limit = float(np.ceil(np.max(np.abs(temperature_panels)) / 2) * 2)


# %% 4. 合并五模式控制相对辐射，固定所有成员计算组均值与组差

join_annual = pd.read_csv(join_source / 'corrected_annual.csv', float_precision='round_trip')
annual = join_annual.loc[~join_annual.model.eq('CAS-ESM2-0')].copy()
cas_annual = pd.read_csv(cas_source / 'annual_global_energy.csv', float_precision='round_trip')
cas_corrected = join_annual.loc[join_annual.model.eq('CAS-ESM2-0')].copy()
cas_control_metadata = json.loads((cas_corrected_source / 'metadata.json').read_text(encoding='utf-8'))
raw_forced_columns = ['model', 'simulation_year', 'phase', 'source_year', 'rsdt', 'rsut', 'rlut', 'asr', 'negative_olr', 'eei']
raw_forced = pd.concat([annual[raw_forced_columns], cas_annual[raw_forced_columns]], ignore_index=True)
assert not raw_forced.duplicated(['model', 'simulation_year']).any()
assert len(raw_forced) == 5 * 280 and raw_forced[['asr', 'negative_olr', 'eei']].notna().all().all()
np.testing.assert_allclose(raw_forced.eei, raw_forced.asr + raw_forced.negative_olr, atol=1e-10)
raw_columns = fields + [name + '_smooth11' for name in fields]
annual = pd.concat([annual, cas_corrected], ignore_index=True)
assert not annual.duplicated(['model', 'simulation_year']).any()
target_index = pd.MultiIndex.from_product([models, np.arange(1, 281)], names=['model', 'simulation_year'])
target = annual.set_index(['model', 'simulation_year']).reindex(target_index)[raw_columns].reset_index()
target['group'] = ''
for group, members in groups.items():
    target.loc[target.model.isin(members), 'group'] = group
    for year in range(1, 281):
        selected = target.loc[target.model.isin(members) & target.simulation_year.eq(year)]
        assert len(selected) == len(members)
        row = {'group': group, 'simulation_year': year, 'target_n_models': len(members),
               'available_n_models': int(selected[fields].notna().all(axis=1).sum())}
        for name in raw_columns:
            # 固定所有成员同时有效才计算，禁止skipna形成两模式恢复组。
            row[name] = float(selected[name].mean()) if selected[name].notna().all() else np.nan
        group_rows.append(row)
group_annual = pd.DataFrame(group_rows)
non_recovery = group_annual.loc[group_annual.group.eq('non_recovery')].set_index('simulation_year')
recovery = group_annual.loc[group_annual.group.eq('recovery')].set_index('simulation_year')
contrast = pd.DataFrame({'simulation_year': np.arange(1, 281)})
for name in raw_columns:
    contrast[name] = non_recovery[name].to_numpy() - recovery[name].to_numpy()
np.testing.assert_allclose(non_recovery.eei_relative_control,
                           non_recovery.asr_relative_control + non_recovery.negative_olr_relative_control,
                           atol=1e-10, rtol=0)
assert group_annual[raw_columns].notna().all().all() and contrast[raw_columns].notna().all().all()
for frame in [group_annual, contrast]:
    for suffix in ['', '_smooth11']:
        np.testing.assert_allclose(frame['eei_relative_control' + suffix],
                                   frame['asr_relative_control' + suffix] + frame['negative_olr_relative_control' + suffix],
                                   atol=1e-10, rtol=0)
radiation_values = np.vstack([group_annual[[name + '_smooth11' for name in fields]].to_numpy(),
                              contrast[[name + '_smooth11' for name in fields]].to_numpy()])
radiation_limits = [float(np.floor(np.nanmin(radiation_values))), float(np.ceil(np.nanmax(radiation_values)))]
for name, frame in [('non_recovery', non_recovery), ('recovery', recovery),
                    ('non_recovery_minus_recovery', contrast.set_index('simulation_year'))]:
    for window, start, stop in [('initial', 1, 10), ('peak_year140', 140, 140), ('late', 271, 280)]:
        for field in fields:
            endpoint_rows.append({'quantity': field, 'series': name, 'window': window,
                                  'value': float(frame.loc[start:stop, field].mean())})


# %% 5. 保存完整组图数据、控制对齐限制与来源

destination.mkdir(parents=True, exist_ok=True)
np.savez_compressed(destination / 'temperature_panels.npz', latitude=latitude, longitude=longitude,
                    temperature_panels=temperature_panels)
pd.DataFrame(membership_rows).to_csv(destination / 'membership.csv', index=False)
target.to_csv(destination / 'model_radiation.csv', index=False)
raw_forced.to_csv(destination / 'model_radiation_raw_complete.csv', index=False)
group_annual.to_csv(destination / 'group_radiation.csv', index=False)
contrast.to_csv(destination / 'radiation_contrast.csv', index=False)
endpoint_table = pd.DataFrame(endpoint_rows)
endpoint_table.loc[endpoint_table.quantity.eq('tas_residual_degC')].to_csv(destination / 'temperature_summary.csv', index=False)
endpoint_table.to_csv(destination / 'endpoint_summary.csv', index=False)
shutil.copy2(temperature_source / 'coastline_global_50m.geojson', destination / 'coastline_global_50m.geojson')
caption = ('Figure3. Top: piControl-relative global net radiation and shortwave/longwave contributions; '
           'bottom: global near-surface air-temperature residuals (years271-280 minus1-10). '
           'Non-recovery contains CESM2 and CNRM-ESM2-1; recovery contains CAS-ESM2-0, GFDL-ESM4 and '
           'NorESM2-LM. Temperature also includes UKESM1-0-LL in the non-recovery group (three versus three models). '
           'Differences are non-recovery minus recovery. All five radiation records use raw forced '
           'fluxes minus fitted piControl backgrounds, without additional initial centering. '
           'The temperature maps show '
           'longitude/latitude labels without internal graticule lines. '
           'No missing values are replaced or partial group averages used. CAS control alignment '
           'uses a provisional zero-counter interpretation; calendar-matched and all-control '
           'fit sensitivities are saved separately because branch calendar origins are inconsistent.')
metadata = {'status': 'complete_with_CAS_alignment_limitation', 'groups': groups, 'group_sizes': [2, 3],
            'temperature_groups': temperature_groups, 'temperature_group_sizes': [3, 3],
            'contrast': 'non_recovery minus recovery', 'radiation_fields': fields,
            'radiation_formula': 'N*=SW*+LW*; SW*=(rsdt-rsut)_forced-control_fit; LW*=-(rlut_forced-control_fit)',
            'radiation_reference': 'Saved piControl-fitted background; no additional initial centering.',
            'radiation_smoothing': 'One centered11year window across the280year sequence, min6only at outer endpoints. GFDL radiation uses actual source1-140 followed by141-280; no overlap or offsets.',
            'radiation_y_limits_Wm2': radiation_limits, 'missing_control_models': [],
            'all_forced_radiation_complete': True,
            'CAS_control_alignment': cas_control_metadata['primary_mapping'],
            'CAS_alignment_limitation': cas_control_metadata['alignment_limitation'],
            'temperature_definition': 'Raw residual2m tas, fixed group members; no piControl correction.',
            'temperature_limit_degC': temperature_limit,
            'temperature_window_alignment': temperature_metadata['window_alignment'],
            'dimensions': 'radiation:model5 xyear280; temperature:panel3 xlat180 xlon360',
            'caption': caption,
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                             for path in source_files]}
(destination / 'metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print(pd.DataFrame(endpoint_rows).to_string(index=False))
print('Five-model corrected radiation complete; CAS alignment sensitivity retained explicitly.')
print(destination.resolve())
