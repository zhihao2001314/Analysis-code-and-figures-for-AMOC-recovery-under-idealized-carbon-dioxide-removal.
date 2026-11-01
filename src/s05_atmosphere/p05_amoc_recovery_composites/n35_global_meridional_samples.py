# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与来源

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
output = base / 'n35_global_meridional_samples'
assert not output.exists(), '正式结果存在，先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']
members = ['r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f2']
starts = [[1, 271], [1, 271], [1, 131], [1850, 2120]]
periods = ['early', 'later']
seasons = ['Annual', 'DJF']
pressure_levels = np.array([1000,2000,3000,5000,7000,10000,15000,20000,25000,30000,40000,50000,60000,70000,85000,92500,100000], dtype=float)
latitude = np.arange(-88., 90., 2.)
radius = 6371000.
gravity = 9.80665
wind_records = json.loads((base / 'n13_endpoint10_moisture_budget/metadata.json').read_text(encoding='utf-8'))['files']
precipitation_records = json.loads((base / 'n01_endpoint_groups/qc.json').read_text(encoding='utf-8'))['inputs']
audit = []
integration_qc = []
reproduction = []
all_models = []
output.mkdir(parents=True)


# %% 3. 每次只读一个完整年的全球月数据，分别构建年平均和DJF样本

for model_index, model in enumerate(models):
    phase_samples = []
    for period_index, period in enumerate(periods):
        first_year = starts[model_index][period_index]
        yearly_samples = []
        for year_index, year in enumerate(range(first_year, first_year + 10)):
            fields = {}
            expected_ids = year * 12 + np.arange(1,13)
            if model != 'NorESM2-LM':
                pieces = []
                paths = sorted(set(row['file'] for row in wind_records if row['model']==model and row['period']==period and row['variable']=='ps'))
                for filename in paths:
                    with xr.open_dataset(filename, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                        take = source.time.dt.year == year
                        if not bool(take.any()):
                            continue
                        assert source.ps.dims == ('time','lat','lon') and source.ps.attrs['units']=='Pa'
                        field = source.ps.isel(time=take).load()
                        ids = field.time.dt.year.values * 12 + field.time.dt.month.values
                        pieces.append(field.assign_coords(time=ids).reset_coords(drop=True))
                        audit.append({'model':model,'period':period,'year':year,'variable':'ps','file':filename,'dimensions':list(source.ps.dims),'months':int(take.sum())})
                surface_pressure = xr.concat(pieces,'time').sortby('time')
                assert np.array_equal(surface_pressure.time.values, expected_ids)

            for variable in ['ua','va','pr']:
                records = precipitation_records if variable=='pr' else wind_records
                paths = sorted(set(row['file'] for row in records if row['model']==model and row['period']==period and row['variable']==variable))
                pieces = []
                embedded_ps = []
                for filename in paths:
                    with xr.open_dataset(filename, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                        take = source.time.dt.year == year
                        if not bool(take.any()):
                            continue
                        assert source.attrs['source_id']==model and source.attrs['variant_label']==members[model_index]
                        if variable=='pr':
                            assert source.pr.dims==('time','lat','lon') and source.pr.attrs['units']=='kg m-2 s-1'
                            field = source.pr.isel(time=take).load().astype(float) * 86400.
                        else:
                            vertical = 'lev' if model=='NorESM2-LM' else 'plev'
                            assert source[variable].dims==('time',vertical,'lat','lon') and source[variable].attrs['units']=='m s-1'
                            native = source[variable].isel(time=take).load().astype(float)
                            if model=='NorESM2-LM':
                                native_ps = source.ps.isel(time=take).load()
                                offset = source.a.values * float(source.p0)
                                fraction = source.b.values
                                converted = np.full((native.sizes['time'],len(pressure_levels),native.sizes['lat'],native.sizes['lon']),np.nan)
                                for month in range(native.sizes['time']):
                                    pressure = offset[:,None,None] + fraction[:,None,None] * native_ps.values[month][None]
                                    wind = native.values[month]
                                    if np.nanmean(pressure[0]) > np.nanmean(pressure[-1]):
                                        pressure = pressure[::-1]
                                        wind = wind[::-1]
                                    assert np.all(np.diff(pressure,axis=0)>0)
                                    for level_index,target in enumerate(pressure_levels):
                                        lower = np.clip(np.sum(pressure<=target,axis=0)-1,0,len(offset)-2)
                                        p0 = np.take_along_axis(pressure,lower[None],axis=0)[0]
                                        p1 = np.take_along_axis(pressure,(lower+1)[None],axis=0)[0]
                                        w0 = np.take_along_axis(wind,lower[None],axis=0)[0]
                                        w1 = np.take_along_axis(wind,(lower+1)[None],axis=0)[0]
                                        valid = (target>=pressure[0]) & (target<=pressure[-1]) & (target<=native_ps.values[month])
                                        converted[month,level_index] = np.where(valid,w0+(w1-w0)*(np.log(target)-np.log(p0))/(np.log(p1)-np.log(p0)),np.nan)
                                field = xr.DataArray(converted,dims=('time','plev','lat','lon'),coords={'time':native.time,'plev':pressure_levels,'lat':native.lat,'lon':native.lon})
                                ids = native.time.dt.year.values*12+native.time.dt.month.values
                                embedded_ps.append(native_ps.assign_coords(time=ids).reset_coords(drop=True))
                            else:
                                field = native.assign_coords(plev=np.rint(native.plev)).sel(plev=pressure_levels)
                        pieces.append(field.reset_coords(drop=True))
                        audit.append({'model':model,'period':period,'year':year,'variable':variable,'file':filename,'dimensions':list(source[variable].dims),'months':int(take.sum())})

                monthly = xr.concat(pieces,'time').sortby('time')
                ids = monthly.time.dt.year.values*12+monthly.time.dt.month.values
                assert np.array_equal(ids,expected_ids)
                month_days = monthly.time.dt.days_in_month.values.astype(float)
                assert np.max(abs(np.diff(monthly.lon.values)-360./monthly.sizes['lon'])) < 1.e-4
                if variable=='ua':
                    reference_days = month_days.copy()
                else:
                    np.testing.assert_array_equal(month_days,reference_days)
                if variable=='pr':
                    assert np.isfinite(monthly.values).all()
                    fields['pr'] = monthly.mean('lon').interp(lat=latitude).assign_coords(time=ids)
                    continue
                if model=='NorESM2-LM':
                    surface_pressure = xr.concat(embedded_ps,'time').sortby('time')
                ps = surface_pressure.sel(time=ids).interp(lat=monthly.lat,lon=monthly.lon).values
                assert np.isfinite(ps).all()
                ps_field = xr.DataArray(ps,dims=('time','lat','lon'),coords={'time':monthly.time,'lat':monthly.lat,'lon':monthly.lon})
                monthly = monthly.where(xr.DataArray(pressure_levels,dims='plev',coords={'plev':pressure_levels})<=ps_field)
                wind = monthly.values
                valid = np.isfinite(wind)
                assert valid[:,0].all() and not np.any((~valid[:,:-1]) & valid[:,1:])
                if variable=='ua':
                    support = xr.DataArray(valid.mean(axis=3)>=.8,dims=('time','plev','lat'),coords={'time':monthly.time,'plev':pressure_levels,'lat':monthly.lat})
                    fields['u'] = monthly.mean('lon',skipna=True).where(support).interp(lat=latitude).assign_coords(time=ids)
                else:
                    filled = wind.copy()
                    for level_index in range(1,len(pressure_levels)):
                        filled[:,level_index] = np.where(valid[:,level_index],filled[:,level_index],filled[:,level_index-1])
                    integral = np.zeros((wind.shape[0],wind.shape[2],wind.shape[3]))
                    profiles = [integral.mean(axis=2).copy()]
                    for level_index in range(len(pressure_levels)-1):
                        width = pressure_levels[level_index+1]-pressure_levels[level_index]
                        thickness = np.clip(ps-pressure_levels[level_index],0.,width)
                        lower_wind = filled[:,level_index]+(filled[:,level_index+1]-filled[:,level_index])*thickness/width
                        integral += .5*(filled[:,level_index]+lower_wind)*thickness
                        profiles.append(integral.mean(axis=2).copy())
                    streamfunction = np.stack(profiles,axis=1) * (2*np.pi*radius*np.cos(np.deg2rad(monthly.lat.values))/gravity)[None,None] / 1.e9
                    field = xr.DataArray(streamfunction,dims=('time','plev','lat'),coords={'time':ids,'plev':pressure_levels,'lat':monthly.lat})
                    fields['psi'] = field.interp(lat=latitude)
                    integration_qc.append({'model':model,'period':period,'year':year,'top_boundary_max_abs':float(abs(streamfunction[:,0]).max()),'maximum_abs_psi_1e9kgs':float(abs(streamfunction).max())})

            monthly_result = xr.Dataset(fields)
            season_fields = []
            for season in seasons:
                take = np.ones(12,dtype=bool) if season=='Annual' else np.isin(np.arange(1,13),[1,2,12])
                weights = xr.DataArray(reference_days[take],dims='time',coords={'time':expected_ids[take]})
                result = monthly_result.isel(time=take).weighted(weights).mean('time',skipna=False)
                result['day_weight'] = float(reference_days[take].sum())
                season_fields.append(result.expand_dims(season=[season]))
            yearly_samples.append(xr.concat(season_fields,'season').expand_dims(year_in_window=[year_index+1]))
            print(model,period,'year',year_index+1,'of10 complete',flush=True)
        phase_samples.append(xr.concat(yearly_samples,'year_in_window').expand_dims(period=[period]))
    result = xr.concat(phase_samples,'period').expand_dims(model=[model])
    result['u'] = result.u.transpose('model','period','season','year_in_window','plev','lat')
    result['psi'] = result.psi.transpose('model','period','season','year_in_window','plev','lat')
    result['pr'] = result.pr.transpose('model','period','season','year_in_window','lat')
    result['day_weight'] = result.day_weight.transpose('model','period','season','year_in_window')
    result.u.attrs['units']='m s-1'
    result.psi.attrs['units']='1e9 kg s-1'
    result.pr.attrs['units']='mm day-1'
    result.attrs={'grouping':'AMOC endpoint recovery, explicitly confirmed by user; not relabeled as deep-convection recovery',
        'psi':'2pi*a*cos(phi)/g times zonal mean column integral of v from10hPa down; terrain-clipped trapezoid; lowest valid v held to surface; no column mass correction',
        'season':'Annual:12months;DJF:Jan/Feb/Dec within calendar year;dayweighted', 'spatial_support':'Native global full longitude; linear interpolation of zonal profiles to88S–88N; no polar extrapolation'}
    with xr.open_dataset(base/'n32_nh_meridional_circulation'/(model+'_samples.nc')) as old:
        for variable in ['u','psi']:
            difference = result[variable].sel(season='DJF',lat=old.lat) - old[variable]
            error = float(abs(difference).max())
            assert error<1.e-8,(model,variable,error)
            reproduction.append({'model':model,'variable':variable,'previous_DJF_section_max_error':error})
    result.to_netcdf(output/(model+'_samples.nc'))
    all_models.append(result)


# %% 4. 保存初始气候态、末期气候态、成员变化及明确分组

samples = xr.concat(all_models,'model')
climate = samples.drop_vars('day_weight').weighted(samples.day_weight).mean('year_in_window',skipna=False)
change = climate.sel(period='later',drop=True)-climate.sel(period='early',drop=True)
groups = xr.concat([climate.sel(model=models[:2]).mean('model',skipna=False).expand_dims(product=['recovered']),
                    climate.sel(model=models[2:]).mean('model',skipna=False).expand_dims(product=['not_recovered'])],'product')
group_change = groups.sel(period='later',drop=True)-groups.sel(period='early',drop=True)
contrast = (group_change.sel(product='not_recovered',drop=True)-group_change.sel(product='recovered',drop=True)).expand_dims(product=['contrast'])
climate.to_netcdf(output/'member_climatology.nc')
change.to_netcdf(output/'member_changes.nc')
groups.to_netcdf(output/'group_climatology.nc')
xr.concat([group_change,contrast],'product').to_netcdf(output/'group_changes.nc')
pd.DataFrame(audit).to_csv(output/'input_audit.csv',index=False)
pd.DataFrame(integration_qc).to_csv(output/'integration_qc.csv',index=False)
pd.DataFrame(reproduction).to_csv(output/'previous_DJF_reproduction.csv',index=False)
(output/'metadata.json').write_text(json.dumps({'models':models,'groups':{'recovered':models[:2],'not_recovered':models[2:]},
    'group_basis':'AMOC late/early>=0.8 after previous initial-state screening; user confirmed reuse, not deep-convection classification',
    'periods':'1pctCO2 1-10 andCDR131-140; no piControl correction', 'seasons':seasons,'starts':starts,
    'dimension_order':{'u_psi':['model','period','season','year_in_window','plev','lat'],'pr':['model','period','season','year_in_window','lat'],'day_weight':['model','period','season','year_in_window']},
    'latitude':latitude.tolist(),'pressure_Pa':pressure_levels.tolist(),'psi':result.attrs['psi'],
    'warning':'Mass circulation diagnostic, not PW energy transport. ITCZ diagnosed from precipitation separately; covariance does not identify causality.'},ensure_ascii=False,indent=2),encoding='utf-8')
print(output,flush=True)
