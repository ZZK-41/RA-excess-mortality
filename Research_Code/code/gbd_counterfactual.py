"""Expected change in GBD RA deaths per person with RA if, in every unit (super-region or country), RA deaths per case had changed in proportion
to that unit's all-cause mortality (constant ratio). Units are combined with their actual case shares in each age-sex cell, so changes in the
regional mix of cases affect observed and expected values equally. Standardised to the global 2019 RA age-sex distribution, ages >=40.
Args: gbd_annual_dir gbd_2019_2023_dir locmap_csv out_dir"""
import sys, io, json, zipfile
from pathlib import Path
import numpy as np, pandas as pd
gbd_annual, gbd_2019_2023, locmap, out = sys.argv[1:5]; out = Path(out)
def rz(p):
    z = zipfile.ZipFile(p); return pd.concat([pd.read_csv(io.BytesIO(z.read(n)), low_memory=False) for n in z.namelist() if n.endswith('.csv')])
AGE40 = [13, 14, 15, 16, 17, 18, 19, 20, 30, 31, 32, 235]; SR = [4, 31, 64, 103, 137, 158, 166]
key = ['location_id', 'sex_id', 'age_id', 'year']
# ---- super-regions, annual
ra = pd.concat([rz(f'{gbd_annual}/IHME-GBD_2023_DATA-{x}-1.zip') for x in ['04c719ed', 'd1117f20', 'f5ad253b']]).drop_duplicates(['measure_id', 'metric_id', 'location_id', 'sex_id', 'age_id', 'year'])
pv = pd.concat([rz(f'{gbd_annual}/IHME-GBD_2023_DATA-d9e5f203-1.zip'), rz(f'{gbd_annual}/IHME-GBD_2023_DATA-a4d3c489-1.zip')]); pv = pv[pv.measure_id == 5]
ac = rz(f'{gbd_annual}/IHME-GBD_2023_DATA-0e1e6620-1.zip')
D = ra[(ra.measure_id == 1) & (ra.metric_id == 1) & ra.sex_id.isin([1, 2]) & ra.age_id.isin(AGE40)].set_index(key).val.rename('d')
C = pv[(pv.metric_id == 1) & pv.sex_id.isin([1, 2]) & pv.age_id.isin(AGE40)].set_index(key).val.rename('c')
M = (ac[ac.sex_id.isin([1, 2]) & ac.age_id.isin(AGE40) & (ac.metric_id == 3)].set_index(key).val/1e5).rename('m')
df = pd.concat([D, C, M], axis=1).dropna().reset_index()
gl = df[df.location_id == 1].set_index(['sex_id', 'age_id', 'year'])
w = gl.xs(2019, level='year').c; w = w/w.sum()
sr = df[df.location_id.isin(SR)].set_index(['location_id', 'sex_id', 'age_id', 'year'])
# global deaths equal the sum of super-regions?
dsum = sr.d.groupby(level=['sex_id', 'age_id', 'year']).sum(); assert ((dsum - gl.d.reindex(dsum.index)).abs()/gl.d.reindex(dsum.index)).max() < 1e-6
a = sr.d/sr.c; base = a.xs(1990, level='year'); mb = sr.m.xs(1990, level='year')
rows = []
for yr in sorted(df.year.unique()):
    s = sr.xs(yr, level='year'); pi = s.c/s.c.groupby(level=['sex_id', 'age_id']).transform('sum')
    obs = (pi*a.xs(yr, level='year')).groupby(level=['sex_id', 'age_id']).sum()
    exp_ = (pi*base*s.m/mb).groupby(level=['sex_id', 'age_id']).sum()
    pi90 = sr.xs(1990, level='year').c; pi90 = pi90/pi90.groupby(level=['sex_id', 'age_id']).transform('sum')
    obs_f = (pi90*a.xs(yr, level='year')).groupby(level=['sex_id', 'age_id']).sum(); exp_f = (pi90*base*s.m/mb).groupby(level=['sex_id', 'age_id']).sum()
    g = gl.xs(yr, level='year')
    rows.append(dict(year=yr, observed=1000*(w*obs).sum(), expected=1000*(w*exp_).sum(), allcause=1000*(w*g.m).sum(),
                     observed_fixed1990=1000*(w*obs_f).sum(), expected_fixed1990=1000*(w*exp_f).sum()))
S = pd.DataFrame(rows).set_index('year')
for c in S.columns: S[c + '_idx'] = 100*S[c]/S.loc[1990, c]
S.to_csv(out/'G10_world_observed_vs_expected_1990_2023.csv')
pc = lambda col, a_, b_: float(100*(S.loc[b_, col]/S.loc[a_, col] - 1))
res = {'superregion_units': {f'{a_}_{b_}': {c: pc(c, a_, b_) for c in ['observed', 'expected', 'allcause', 'observed_fixed1990', 'expected_fixed1990']}
                             for a_, b_ in [(1990, 2023), (1990, 2019), (2000, 2023)]}}
for k, v in res['superregion_units'].items():
    v['share_of_observed_decline_expected_log'] = float(np.log(1 + v['expected']/100)/np.log(1 + v['observed']/100))
    v['observed_over_expected_pct'] = float(100*((1 + v['observed']/100)/(1 + v['expected']/100) - 1))
# ---- countries, 2000 -> 2023 (also 2010 -> 2023)
R00, A00, R10 = rz(f'{gbd_annual}/IHME-GBD_2023_DATA-80a09abe-1.zip'), rz(f'{gbd_annual}/IHME-GBD_2023_DATA-532d5fea-1.zip'), rz(f'{gbd_annual}/IHME-GBD_2023_DATA-4af37d63-1.zip')
R23, A1923 = rz(f'{gbd_2019_2023}/IHME-GBD_2023_DATA-aaf6002d-1.zip'), rz(f'{gbd_2019_2023}/IHME-GBD_2023_DATA-8e54c028-1.zip')
lm = pd.read_csv(locmap, encoding='utf-8-sig'); countries = set(lm.location_id)
k3 = ['location_id', 'sex_id', 'age_id']
def cells(r, acx, yr):
    x = r[(r.metric_id == 1) & r.sex_id.isin([1, 2]) & r.age_id.isin(AGE40) & (r.year == yr) & r.location_id.isin(countries)].pivot_table(index=k3, columns='measure_id', values='val').rename(columns={1: 'd', 5: 'c'})
    m = acx[acx.sex_id.isin([1, 2]) & acx.age_id.isin(AGE40) & (acx.year == yr) & (acx.metric_id == 3) & acx.location_id.isin(countries)].set_index(k3).val.rename('m')/1e5
    return x[['d', 'c']].join(m, how='inner')
c00, c10, c23 = cells(R00, A00, 2000), cells(R10, A00, 2010), cells(R23, A1923, 2023)
assert c00.index.get_level_values(0).nunique() == 204 and c23.index.get_level_values(0).nunique() == 204
def std(cc, rate):  # case-weighted across countries within cell, then standardised
    pi = cc.c/cc.c.groupby(level=['sex_id', 'age_id']).transform('sum')
    return 1000*(w*(pi*rate).groupby(level=['sex_id', 'age_id']).sum()).sum()
cn = {}
for lab, c0 in [('2000_2023', c00), ('2010_2023', c10)]:
    o0 = std(c0, c0.d/c0.c); o1 = std(c23, c23.d/c23.c)
    e1 = std(c23, (c0.d/c0.c).reindex(c23.index)*c23.m/c0.m.reindex(c23.index))
    cn[lab] = dict(observed=float(100*(o1/o0 - 1)), expected=float(100*(e1/o0 - 1)))
    cn[lab]['share_of_observed_decline_expected_log'] = float(np.log(1 + cn[lab]['expected']/100)/np.log(1 + cn[lab]['observed']/100))
    cn[lab]['observed_over_expected_pct'] = float(100*((1 + cn[lab]['observed']/100)/(1 + cn[lab]['expected']/100) - 1))
res['country_units'] = cn
json.dump(res, open(out/'G10_counterfactual.json', 'w'), indent=1)
print(json.dumps(res, indent=1, default=lambda x: round(x, 2)))
print(S[['observed', 'expected', 'allcause', 'observed_idx', 'expected_idx', 'allcause_idx']].iloc[::3].round(2).to_string())
