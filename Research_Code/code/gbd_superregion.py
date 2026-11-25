"""GBD 2023, ages 45-89: deaths attributed to RA vs modelled excess deaths under the four mortality profiles,
by GBD super-region. Args: location_map_csv cells_csv_gz out_csv"""
import sys
import numpy as np, pandas as pd
locmap, cells_path, out = sys.argv[1:4]
mapping = pd.read_csv(locmap, encoding='utf-8-sig').set_index('location_id').super_region
cl = pd.read_csv(cells_path, low_memory=False)
cl = cl[(cl.location_name != 'Global') & (cl.scenario != 'Canada 2019 fixed') & (cl.year == 2023)].copy()
cl['super_region'] = cl.location_id.map(mapping)
print('unmapped locations:', cl[cl.super_region.isna()].location_name.unique()[:10], cl.location_id.nunique())
base = cl[cl.scenario == 'Canada relative'].groupby('super_region').agg(cases=('cases', 'sum'), gbd_coded=('ra_deaths', 'sum'))
ex = cl.pivot_table(index='super_region', columns='scenario', values='equivalent', aggfunc='sum')
r = base.join(ex)
r['coded_per_1000_cases'] = 1000*r.gbd_coded/r.cases
for s in ex.columns: r[f'ratio_{s}'] = r[s]/r.gbd_coded
r['ratio_min'] = r[[f'ratio_{s}' for s in ex.columns]].min(axis=1); r['ratio_max'] = r[[f'ratio_{s}' for s in ex.columns]].max(axis=1)
tot = r[['cases', 'gbd_coded'] + list(ex.columns)].sum(); tot['coded_per_1000_cases'] = 1000*tot.gbd_coded/tot.cases
for s in ex.columns: tot[f'ratio_{s}'] = tot[s]/tot.gbd_coded
tot['ratio_min'] = min(tot[f'ratio_{s}'] for s in ex.columns); tot['ratio_max'] = max(tot[f'ratio_{s}'] for s in ex.columns)
r.loc['All 204 locations'] = tot
r = r.sort_values('coded_per_1000_cases'); r.to_csv(out)
pd.set_option('display.width', 250)
print(r[['cases', 'gbd_coded', 'coded_per_1000_cases', 'ratio_min', 'ratio_max']].round(2).to_string())
# country-level coded deaths per 1000 cases: spread
c = cl[cl.scenario == 'Canada relative'].groupby('location_name').agg(cases=('cases', 'sum'), coded=('ra_deaths', 'sum'))
c['per1000'] = 1000*c.coded/c.cases
print('country coded per 1000 cases: median', round(c.per1000.median(), 2), 'IQR', c.per1000.quantile([.25, .75]).round(2).tolist(), 'min-max', round(c.per1000.min(), 2), round(c.per1000.max(), 2))
print(c.loc[['Canada', 'France', 'Denmark', 'Italy', 'United States of America', 'United Kingdom', 'China', 'India', 'Japan', 'Brazil']].round(2).to_string())
