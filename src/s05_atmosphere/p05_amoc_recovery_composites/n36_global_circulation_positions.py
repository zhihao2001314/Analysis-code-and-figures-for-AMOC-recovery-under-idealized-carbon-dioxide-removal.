# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.integrate import cumulative_trapezoid, trapezoid


# %% 2. 路径与操作性指标定义

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
source_directory = base / 'n35_global_meridional_samples'
output = base / 'n36_global_circulation_positions'
assert not output.exists(), '正式结果须先归档。'
models = ['GFDL-ESM4','NorESM2-LM','CESM2','CNRM-ESM2-1']
seasons = ['Annual','DJF']
itcz_bounds = [-20.,20.]
hadley_peak_bound = 30.
hadley_edge_search_limit = 55.
jet_definitions = [('subtropical_u250',25000,15.,45.),('midlatitude_u850',85000,30.,70.)]
rows = []


# %% 3. 读取成员和组气候态，不用成员位置的平均冒充组位置

with xr.open_dataset(source_directory/'member_climatology.nc') as source:
    member_climate = source.load()
with xr.open_dataset(source_directory/'group_climatology.nc') as source:
    group_climate = source.load()
states = xr.concat([member_climate.rename(model='product'),group_climate],'product')
latitude = states.lat.values


# %% 4. 两种ITCZ定位、Hadley边界和高低层西风峰值

for product in states['product'].values:
    for season in seasons:
        for period in ['early','later']:
            state = states.sel(product=product,season=season,period=period)
            tropical = state.pr.sel(lat=slice(*itcz_bounds))
            lat = tropical.lat.values
            rainfall = tropical.values
            weights = rainfall*np.cos(np.deg2rad(lat))
            cumulative = cumulative_trapezoid(weights,lat,initial=0.)
            total = cumulative[-1]
            assert np.isfinite(weights).all() and total>0
            centroid = float(trapezoid(lat*weights,lat)/total)
            median = float(np.interp(total/2,cumulative,lat))
            peak = float(lat[np.argmax(rainfall)])
            for metric,value in [('ITCZ_precipitation_median',median),('ITCZ_precipitation_centroid',centroid),('tropical_precipitation_peak_grid',peak)]:
                rows.append({'product':str(product),'season':season,'period':period,'metric':metric,'value':value,'unit':'degree_north','status':'diagnosed'})

            psi = state.psi.sel(plev=50000)
            rows.append({'product':str(product),'season':season,'period':period,'metric':'cross_equatorial_psi500','value':float(psi.sel(lat=0)),'unit':'1e9 kg s-1','status':'diagnosed'})
            for hemisphere,sign in [('NH',1.),('SH',-1.)]:
                poleward_latitude = latitude[(sign*latitude>=0)&(sign*latitude<=hadley_edge_search_limit)]
                poleward_latitude = poleward_latitude[np.argsort(sign*poleward_latitude)]
                profile = psi.sel(lat=poleward_latitude).values*sign
                tropical_indices = np.flatnonzero(sign*poleward_latitude<=hadley_peak_bound)
                maximum_index = int(tropical_indices[np.argmax(profile[tropical_indices])])
                edge = np.nan
                status = 'no_positive_direct_cell_or_no_poleward_zero_crossing'
                if profile[maximum_index]>0:
                    crossings = np.flatnonzero((np.arange(len(profile)-1)>=maximum_index)&(profile[:-1]>0)&(profile[1:]<=0))
                    if len(crossings)>0:
                        index = int(crossings[0])
                        edge = float(poleward_latitude[index]-profile[index]*(poleward_latitude[index+1]-poleward_latitude[index])/(profile[index+1]-profile[index]))
                        status = 'diagnosed'
                rows.append({'product':str(product),'season':season,'period':period,'metric':hemisphere+'_Hadley_edge_psi500','value':edge,'unit':'degree_north','status':status})
                strength = float(psi.sel(lat=poleward_latitude[maximum_index]))
                rows.append({'product':str(product),'season':season,'period':period,'metric':hemisphere+'_Hadley_psi500_extremum','value':strength,'unit':'1e9 kg s-1','status':'diagnosed' if profile[maximum_index]>0 else 'direct_cell_not_identified'})

                for name,pressure,equatorward,poleward in jet_definitions:
                    selected = state.u.sel(plev=pressure).where((sign*state.lat>=equatorward)&(sign*state.lat<=poleward),drop=True)
                    values = selected.values
                    lat = selected.lat.values
                    position = np.nan
                    speed = np.nan
                    status = 'no_interior_positive_peak'
                    if np.isfinite(values).any():
                        index = int(np.nanargmax(values))
                        if 0<index<len(values)-1 and values[index]>0 and np.isfinite(values[index-1:index+2]).all():
                            curvature = values[index-1]-2*values[index]+values[index+1]
                            shift = .5*(values[index-1]-values[index+1])/curvature
                            assert abs(shift)<=1.
                            position = float(lat[index]+shift*(lat[index+1]-lat[index]))
                            speed = float(values[index]-.25*(values[index-1]-values[index+1])*shift)
                            status = 'three_point_parabolic_peak_on_2degree_profile'
                    rows.append({'product':str(product),'season':season,'period':period,'metric':hemisphere+'_'+name+'_latitude','value':position,'unit':'degree_north','status':status})
                    rows.append({'product':str(product),'season':season,'period':period,'metric':hemisphere+'_'+name+'_speed','value':speed,'unit':'m s-1','status':status})


# %% 5. 保存位置与端点变化，缺失的边界峰值不强行解释

table = pd.DataFrame(rows)
changes = table.pivot(index=['product','season','metric','unit'],columns='period',values='value').reset_index()
changes['later_minus_early'] = changes['later']-changes['early']
output.mkdir(parents=True)
table.to_csv(output/'endpoint_positions.csv',index=False)
changes.to_csv(output/'position_changes.csv',index=False)
metadata = {'groups':{'recovered':models[:2],'not_recovered':models[2:]},'group_basis':'AMOC, confirmed by user',
    'ITCZ':'Precipitation area-weighted median and first-moment centroid within20S–20N, all longitudes and all surfaces; grid maximum only as sensitivity. Neither is a direct diagnosed vertical velocity.',
    'Hadley_edge':'500hPa streamfunction first poleward zero crossing after sign-consistent tropical extremum; no forced edge when direct cell absent.',
    'jets':jet_definitions,'jet_method':'Positive interior maximum among valid latitudes, three finite neighboring points for parabolic refinement of2degree zonal profile. Some low-level latitudes fail80percent longitude support; no filling across those gaps. Operational upper/lower wind maxima, not unique dynamically separated jet axes.',
    'sign':'Latitude is degrees north; positive later-minus-early means northward for both hemispheres, not poleward in both.',
    'nonlinearity':'Metrics diagnosed from each group-mean climatology separately, not averaged member metrics.',
    'inference':'Exploratory endpoint positions, no significance test of shifts in this first-look product; no causal identification.',
    'dimensions':'table rows:product,season,period,metric; values and units explicit'}
(output/'metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
print(changes[(changes['product'].isin(['recovered','not_recovered']))&changes.metric.isin(['ITCZ_precipitation_median','NH_Hadley_edge_psi500','NH_subtropical_u250_latitude','NH_midlatitude_u850_latitude'])].round(4).to_string(index=False),flush=True)
