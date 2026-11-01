# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
output = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n13_endpoint10_moisture_budget'
assert not output.exists(), '已有正式输出须先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']
members = {'GFDL-ESM4': 'r1i1p1f1', 'NorESM2-LM': 'r1i1p1f1', 'CESM2': 'r1i1p1f1', 'CNRM-ESM2-1': 'r1i1p1f2'}
starts = {'GFDL-ESM4': [1, 271], 'NorESM2-LM': [1, 271], 'CESM2': [1, 131], 'CNRM-ESM2-1': [1850, 2120]}
gravity = 9.80665
radius = 6371000.
top_pressure = 1000.
periods = ['early', 'later']
nor_pr = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n05_nor_pr_monthly_from_daily/nor_pr_cdr_027101-028012_monthly_derived.nc'
audit = []
integration_qc = []
monthly_rows = []
output.mkdir(parents=True)


# %% 3. 逐模式读取整窗及时间导数邻月，科学口径为30个DJF日历月

for model in models:
    phases = []
    phase_maps = []
    for phase_index, period in enumerate(periods):
        experiment = '1pctCO2' if period == 'early' else '1pctCO2-cdr'
        start = starts[model][phase_index]
        lower_id = start * 12
        upper_id = (start + 10) * 12 + 1
        fields = {}
        for variable in ['hus', 'ua', 'va']:
            table = ('CFmon' if variable == 'hus' else 'AERmon') if model == 'NorESM2-LM' else 'Amon'
            folder = Path('F:/sea/cmip6' + variable + ('ad' if period == 'early' else ''))
            paths = sorted(folder.rglob(f'{variable}_{table}_{model}_{experiment}_{members[model]}_*.nc'))
            paths = [path for path in paths if '(' not in path.name and 'yearmean' not in path.name.lower()]
            parts = []
            ps_parts = []
            for path in paths:
                date_range = path.stem.split('_')[-1]
                if int(date_range.split('-')[0][:4]) > start + 10 or int(date_range.split('-')[1][:4]) < start - 1:
                    continue
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    ids = source.time.dt.year * 12 + source.time.dt.month
                    take = (ids >= lower_id) & (ids <= upper_id)
                    if not bool(take.any()):
                        continue
                    assert source.attrs['source_id'] == model and source.attrs['experiment_id'] == experiment
                    assert source.attrs['variant_label'] == members[model]
                    vertical = 'lev' if model == 'NorESM2-LM' else 'plev'
                    assert source[variable].dims == ('time', vertical, 'lat', 'lon')
                    assert source[variable].attrs['units'] == ('1' if variable == 'hus' else 'm s-1')
                    print('Reading', path.name, int(take.sum()), 'months', flush=True)
                    # 先按原生经度顺序读连续纬带，再重排经度，避免HDF压缩块被逐列重复解压。
                    latitude_take = (source.lat >= 15) & (source.lat <= 85)
                    field = source[variable].isel(time=take, lat=latitude_take).load()
                    field = field.assign_coords(lon=(source.lon + 180) % 360 - 180).sortby('lon').sortby('lat').sel(lon=slice(-90, 40)).astype('float64')
                    if vertical == 'plev':
                        field = field.assign_coords(plev=np.rint(field.plev)).sortby('plev').sel(plev=slice(top_pressure, 100000.))
                        assert float(field.plev[0]) == top_pressure
                    else:
                        assert source.ps.dims == ('time', 'lat', 'lon') and source.ps.attrs['units'] == 'Pa'
                        native_ps = source.ps.isel(time=take, lat=latitude_take).load().assign_coords(lon=(source.lon + 180) % 360 - 180).sortby('lon').sortby('lat').sel(lat=field.lat, lon=field.lon).astype('float64')
                        ps_parts.append(native_ps)
                        if variable == 'hus':
                            native_coefficients = {name: source[name].values.copy() for name in ['a_bnds', 'b_bnds', 'a', 'b', 'lev', 'p0']}
                        else:
                            for name, value in native_coefficients.items():
                                np.testing.assert_allclose(source[name].values, value, atol=1e-10, rtol=0)
                    parts.append(field.reset_coords(drop=True))
                    audit.append({'model': model, 'period': period, 'variable': variable, 'file': str(path), 'dimensions': list(source[variable].dims), 'units': source[variable].attrs['units'], 'calendar': source.time.encoding.get('calendar', 'standard'), 'selected_months_with_halo': int(take.sum())})
            assert parts, (model, period, variable)
            monthly = xr.concat(parts, dim='time', join='exact').sortby('time')
            ids = monthly.time.dt.year.values * 12 + monthly.time.dt.month.values
            assert np.all(np.diff(ids) == 1) and len(ids) in [121, 122], (model, period, variable, ids)
            if variable == 'hus':
                reference_ids = ids
                dates = monthly.time.values
                days = monthly.time.dt.days_in_month.values.astype(float)
                months = monthly.time.dt.month.values
                target_lat = monthly.lat
                target_lon = monthly.lon
            else:
                assert np.array_equal(ids, reference_ids)
            fields[variable] = monthly.assign_coords(time=ids)
            if model == 'NorESM2-LM':
                native_ps = xr.concat(ps_parts, dim='time', join='exact').sortby('time').assign_coords(time=ids)
                if variable == 'hus':
                    fields['ps'] = native_ps
                else:
                    np.testing.assert_allclose(native_ps.values, fields['ps'].values, atol=.1, rtol=0)
            print(model, period, variable, len(ids), 'months including derivative halo', flush=True)
        if model != 'NorESM2-LM':
            parts = []
            folder = Path('F:/sea/cmip6ps' + ('ad' if period == 'early' else ''))
            for path in sorted(folder.rglob(f'ps_Amon_{model}_{experiment}_{members[model]}_*.nc')):
                if '(' in path.name or 'yearmean' in path.name.lower():
                    continue
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    ids = source.time.dt.year * 12 + source.time.dt.month
                    take = (ids >= lower_id) & (ids <= upper_id)
                    if not bool(take.any()):
                        continue
                    assert source.ps.dims == ('time', 'lat', 'lon') and source.ps.attrs['units'] == 'Pa'
                    field = source.ps.isel(time=take).load().assign_coords(lon=(source.lon + 180) % 360 - 180).sortby('lon').sortby('lat').astype('float64')
                    field = field.interp(lat=target_lat, lon=target_lon)
                    parts.append(field.reset_coords(drop=True))
                    audit.append({'model': model, 'period': period, 'variable': 'ps', 'file': str(path), 'dimensions': list(source.ps.dims), 'selected_months_with_halo': int(take.sum())})
            ps = xr.concat(parts, dim='time', join='exact').sortby('time')
            ids = ps.time.dt.year.values * 12 + ps.time.dt.month.values
            assert np.array_equal(ids, reference_ids)
            fields['ps'] = ps.assign_coords(time=ids)
        q, u, v, ps = xr.align(fields['hus'], fields['ua'], fields['va'], fields['ps'], join='exact')
        q, u, v, ps = q.values.copy(), u.values.copy(), v.values.copy(), ps.values
        assert np.isfinite(ps).all()


        # %% 4. 统一10hPa至地面积分，比湿为水汽而非云水

        if model == 'NorESM2-LM':
            assert np.isfinite(q).all() and np.isfinite(u).all() and np.isfinite(v).all()
            bounds = native_coefficients['a_bnds'][None, :, :, None, None] * float(native_coefficients['p0']) + native_coefficients['b_bnds'][None, :, :, None, None] * ps[:, None, None]
            bottom = np.max(bounds, axis=2)
            top = np.min(bounds, axis=2)
            thickness = np.maximum(bottom - np.maximum(top, top_pressure), 0.)
            mass = thickness / gravity
        else:
            pressure = fields['hus'].plev.values
            valid = np.isfinite(q) & np.isfinite(u) & np.isfinite(v) & (pressure[None, :, None, None] <= ps[:, None])
            assert valid[:, 0].all()
            assert not np.any((~valid[:, :-1]) & valid[:, 1:]), '垂直缺测不连续，停止计算。'
            mass = np.zeros_like(q)
            for k in range(len(pressure) - 1):
                thickness = (pressure[k + 1] - pressure[k]) * valid[:, k + 1]
                mass[:, k] += thickness / 2
                mass[:, k + 1] += thickness / 2
            for k in range(len(pressure)):
                last = valid[:, k] & (~valid[:, k + 1] if k + 1 < len(pressure) else True)
                mass[:, k] += np.where(last, ps - pressure[k], 0.)
                if k > 0:
                    # 地下延拓仅定义分解基准；月积分对应层权重为零。
                    q[:, k] = np.where(valid[:, k], q[:, k], q[:, k - 1])
                    u[:, k] = np.where(valid[:, k], u[:, k], u[:, k - 1])
                    v[:, k] = np.where(valid[:, k], v[:, k], v[:, k - 1])
            mass /= gravity
        mass_error = float(np.max(np.abs(mass.sum(axis=1) * gravity - (ps - top_pressure))))
        assert mass_error < 1e-6
        assert np.isfinite(q).all() and np.isfinite(u).all() and np.isfinite(v).all()
        water = np.sum(q * mass, axis=1)
        flux_u = np.sum(q * u * mass, axis=1)
        flux_v = np.sum(q * v * mass, axis=1)
        seconds = np.array([(date - dates[0]).total_seconds() for date in dates])
        tendency = np.gradient(water, seconds, axis=0, edge_order=2) * 86400.
        inside = (reference_ids >= start * 12 + 1) & (reference_ids <= (start + 10) * 12)
        selected = inside & np.isin(months, [12, 1, 2])
        assert int(inside.sum()) == 120 and int(selected.sum()) == 30
        time_weights = np.zeros(len(reference_ids))
        time_weights[selected] = days[selected] / days[selected].sum()
        w3 = time_weights[:, None, None]
        w4 = time_weights[:, None, None, None]
        state = {'q': np.sum(q*w4, axis=0), 'u': np.sum(u*w4, axis=0), 'v': np.sum(v*w4, axis=0), 'mass': np.sum(mass*w4, axis=0), 'flux_u': np.sum(flux_u*w3, axis=0), 'flux_v': np.sum(flux_v*w3, axis=0), 'storage': np.sum(tendency*w3, axis=0), 'water': np.sum(water*w3, axis=0)}
        column = xr.Dataset({'water': (('time', 'lat', 'lon'), water), 'flux_u': (('time', 'lat', 'lon'), flux_u), 'flux_v': (('time', 'lat', 'lon'), flux_v), 'storage': (('time', 'lat', 'lon'), tendency), 'djf_weight': ('time', time_weights)}, coords={'time': dates, 'lat': target_lat, 'lon': target_lon})
        column.attrs.update({'model': model, 'period': period, 'top_pressure_Pa': top_pressure, 'storage_unit': 'mm day-1', 'water_unit': 'kg m-2', 'flux_unit': 'kg m-1 s-1'})
        column.to_netcdf(output / f'{model}_{period}_columns.nc')
        integration_qc.append({'model': model, 'period': period, 'mass_pressure_max_error_Pa': mass_error, 'selected_DJF_months': 30, 'months_with_halo': len(reference_ids), 'storage_boundary': 'One-sided three-point at early first January; other selected months have neighboring monthly centers.', 'native_integration': 'Hybrid layer dp/g, clipped at 10hPa' if model == 'NorESM2-LM' else 'Pressure-level trapezoid, lowest valid level held to surface, 10hPa top'})


        # %% 5. 同窗降水蒸发重新读取，与既有组差口径对照

        for variable in ['pr', 'evspsbl']:
            folder = Path('F:/sea/cmip6' + variable + ('ad' if period == 'early' else ''))
            paths = sorted(folder.rglob(f'{variable}_Amon_{model}_{experiment}_{members[model]}_*.nc'))
            paths = [path for path in paths if '(' not in path.name and 'yearmean' not in path.name.lower()]
            if model == 'NorESM2-LM' and period == 'later' and variable == 'pr':
                paths = [nor_pr]
            parts = []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    ids = source.time.dt.year * 12 + source.time.dt.month
                    take = (ids >= start * 12 + 1) & (ids <= (start + 10) * 12)
                    if not bool(take.any()):
                        continue
                    assert source[variable].dims == ('time', 'lat', 'lon')
                    assert source[variable].attrs['units'] == 'kg m-2 s-1'
                    field = source[variable].isel(time=take).load().assign_coords(lon=(source.lon + 180) % 360 - 180).sortby('lon').sortby('lat').astype('float64')
                    field = field.interp(lat=target_lat, lon=target_lon)
                    parts.append(field.reset_coords(drop=True))
                    audit.append({'model': model, 'period': period, 'variable': variable, 'file': str(path), 'dimensions': list(source[variable].dims), 'selected_months': int(take.sum())})
            surface = xr.concat(parts, dim='time', join='exact').sortby('time')
            surface_ids = surface.time.dt.year.values * 12 + surface.time.dt.month.values
            assert np.array_equal(surface_ids, reference_ids[inside])
            state[variable] = np.sum(surface.values * time_weights[inside, None, None], axis=0) * 86400.
        phases.append(state)
        print(model, period, 'column transport and surface budget calculated', flush=True)


    # %% 6. 分解两个端点的输送，原生网格计算球面辐合和残差

    latitude, longitude = target_lat.values, target_lon.values
    phi, lam = np.deg2rad(latitude), np.deg2rad(longitude)
    cosine = np.cos(phi)[:, None]
    early, later = phases
    dq = later['q'] - early['q']
    dm = later['mass'] - early['mass']
    flux_terms = {}
    for component in ['u', 'v']:
        dv = later[component] - early[component]
        terms = {'total': later['flux_' + component] - early['flux_' + component], 'dynamic': np.sum(early['q'] * dv * early['mass'], axis=0), 'humidity': np.sum(dq * early[component] * early['mass'], axis=0), 'nonlinear': np.sum(dq * dv * early['mass'], axis=0), 'pressure': np.sum(later['q'] * later[component] * dm, axis=0)}
        terms['monthly_covariance'] = terms['total'] - terms['dynamic'] - terms['humidity'] - terms['nonlinear'] - terms['pressure']
        flux_terms[component] = terms
    result = xr.Dataset(coords={'lat': latitude, 'lon': longitude})
    for term in flux_terms['u']:
        east, north = flux_terms['u'][term], flux_terms['v'][term]
        convergence = -(np.gradient(east, lam, axis=1, edge_order=2) + np.gradient(north * cosine, phi, axis=0, edge_order=2)) / (radius * cosine) * 86400.
        result['mfc_' + term] = (('lat', 'lon'), convergence, {'units': 'mm day-1'})
        for component in ['u', 'v']:
            result[f'flux_{term}_{component}'] = (('lat', 'lon'), flux_terms[component][term], {'units': 'kg m-1 s-1'})
    for name in ['pr', 'evspsbl', 'storage', 'water']:
        result['delta_' + name] = (('lat', 'lon'), later[name] - early[name])
    result['delta_pe'] = result.delta_pr - result.delta_evspsbl
    result['resolved_pe'] = result.mfc_total - result.delta_storage
    result['budget_residual'] = result.delta_pe - result.resolved_pe
    result['mfc_other'] = result.mfc_nonlinear + result.mfc_pressure + result.mfc_monthly_covariance
    result['flux_other_u'] = result.flux_nonlinear_u + result.flux_pressure_u + result.flux_monthly_covariance_u
    result['flux_other_v'] = result.flux_nonlinear_v + result.flux_pressure_v + result.flux_monthly_covariance_v
    decomposition_error = float(abs(result.mfc_total - result.mfc_dynamic - result.mfc_humidity - result.mfc_other).max())
    assert decomposition_error < 1e-9
    climate = []
    for period, state in zip(periods, phases):
        maps = xr.Dataset({name: (('lat', 'lon'), state[name]) for name in ['pr', 'evspsbl', 'storage', 'water', 'flux_u', 'flux_v']}, coords={'lat': latitude, 'lon': longitude})
        climate.append(maps.expand_dims(period=[period]))
    xr.concat(climate, dim='period').sel(lat=slice(20, 80), lon=slice(-80, 30)).to_netcdf(output / f'{model}_climatology.nc')
    result = result.sel(lat=slice(20, 80), lon=slice(-80, 30))
    result.attrs.update({'model': model, 'endpoint': 'CDR131-140 minus 1pctCO2 1-10; Jan/Feb/Dec within each ten-calendar-year window, 30 months, day weighted', 'top_pressure_Pa': top_pressure, 'integration': 'Nor native hybrid layer masses; other models pressure-level trapezoid. Same 10hPa top; differing vertical representations retained.', 'decomposition': 'Early-reference dynamic + humidity + nonlinear + pressure + monthly covariance; monthly covariance is not submonthly transient eddy flux.', 'residual': 'delta(P-E) - (monthly mean-flow MFC change - column-water storage change). Includes unresolved submonthly transport, sampling and numerical/vertical integration errors.', 'decomposition_error_mm_day': decomposition_error})
    result.to_netcdf(output / f'{model}.nc')
    print(model, 'saved, decomposition max error', decomposition_error, flush=True)


# %% 7. 保存输入来源与统一口径

metadata = {'models': models, 'groups': {'recovered': models[:2], 'not_recovered': models[2:]}, 'files': audit, 'integration_qc': integration_qc, 'window': 'later CDR131-140 minus early1pctCO2 1-10; 30 DJF calendar months per period; day weighted; not nine complete winters', 'dimensions': {'raw_pressure_fields': ['time', 'plev', 'lat', 'lon'], 'raw_Nor_fields': ['time', 'lev', 'lat', 'lon'], 'column_products': ['time', 'lat', 'lon'], 'native_changes': ['lat', 'lon'], 'climatology': ['period', 'lat', 'lon']}, 'top_pressure_Pa': top_pressure, 'gravity': gravity, 'earth_radius_m': radius, 'drift_correction': False, 'surface_data': 'Nor CDR precipitation reuses verified daily-to-monthly derived series. Others published monthly.', 'limitations': 'Monthly products omit submonthly covariance; no closed full moisture-budget or causal attribution claimed. Nor native-layer versus other pressure-level integration is a methodological difference, not automatically an improvement in closure.'}
(output / 'metadata.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding='utf-8')
print('OUTPUT', output, flush=True)
