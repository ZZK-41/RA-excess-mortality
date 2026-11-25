"""Observed GBD 2023 RA deaths per person with RA (ages >=40, standardised to the global 2019 RA age-sex distribution) and the values
expected if, in every location and age-sex group, RA deaths per person with RA had changed in proportion to all-cause mortality
(constant relative excess mortality and constant certification). Locations are combined with their actual numbers of cases in each year,
so changes in the mix of locations affect observed and expected values equally.
  World annual 1990-2023: units = 7 super-regions.  World and super-regions 2000-2023: units = 204 countries and territories.
  Super-regions 1990-2023: the super-region is one unit.  Countries: the country is the unit.
Args: gbd_annual_dir gbd_2019_2023_dir locmap_csv out_dir figdata_dir"""
import sys, io, json, zipfile
from pathlib import Path
import numpy as np, pandas as pd
gbd_annual, gbd_2019_2023, locmap, out, fd = sys.argv[1:6]; out, fd = Path(out), Path(fd)
def rz(p):
    z = zipfile.ZipFile(p); return pd.concat([pd.read_csv(io.BytesIO(z.read(n)), low_memory=False) for n in z.namelist() if n.endswith('.csv')])
AGE40 = [13, 14, 15, 16, 17, 18, 19, 20, 30, 31, 32, 235]; SR = [4, 31, 64, 103, 137, 158, 166]
NAMES = {1: 'Global', 64: 'High-income', 31: 'Central Europe, Eastern Europe, and Central Asia', 103: 'Latin America and Caribbean', 137: 'North Africa and Middle East',
         158: 'South Asia', 4: 'Southeast Asia, East Asia, and Oceania', 166: 'Sub-Saharan Africa', 101: 'Canada', 78: 'Denmark', 80: 'France', 86: 'Italy',
         71: 'Australia', 79: 'Finland', 67: 'Japan', 89: 'Netherlands', 93: 'Sweden', 95: 'United Kingdom', 102: 'United States of America'}
HI11 = [71, 101, 78, 79, 80, 86, 67, 89, 93, 95, 102]
key = ['location_id', 'sex_id', 'age_id', 'year']; cell = ['sex_id', 'age_id']
# ---------------- annual data, 19 locations
ra = pd.concat([rz(f'{gbd_annual}/IHME-GBD_2023_DATA-{x}-1.zip') for x in ['04c719ed', 'd1117f20', 'f5ad253b']]).drop_duplicates(['measure_id', 'metric_id', 'location_id', 'sex_id', 'age_id', 'year'])
pv = pd.concat([rz(f'{gbd_annual}/IHME-GBD_2023_DATA-d9e5f203-1.zip'), rz(f'{gbd_annual}/IHME-GBD_2023_DATA-a4d3c489-1.zip')]); pv = pv[pv.measure_id == 5]
ac = rz(f'{gbd_annual}/IHME-GBD_2023_DATA-0e1e6620-1.zip')
D = ra[(ra.measure_id == 1) & (ra.metric_id == 1) & ra.sex_id.isin([1, 2]) & ra.age_id.isin(AGE40)].set_index(key).val.rename('d')
C = pv[(pv.metric_id == 1) & pv.sex_id.isin([1, 2]) & pv.age_id.isin(AGE40)].set_index(key).val.rename('c')
M = (ac[ac.sex_id.isin([1, 2]) & ac.age_id.isin(AGE40) & (ac.metric_id == 3)].set_index(key).val/1e5).rename('m')
an = pd.concat([D, C, M], axis=1).dropna()
assert an.groupby(level=['location_id', 'year']).size().eq(24).all()
w = an.xs((1, 2019), level=['location_id', 'year']).c; w = w/w.sum()
def series(units, t0, data):
    """units: list of location ids combined by case shares; returns DataFrame year -> observed, expected (per 1000, standardised)."""
    x = data[data.index.get_level_values('location_id').isin(units)]
    base = x.xs(t0, level='year'); a0 = base.d/base.c; m0 = base.m
    rows = []
    for yr, g in x.groupby(level='year'):
        g = g.droplevel('year'); pi = g.c/g.c.groupby(level=cell).transform('sum')
        obs = (pi*g.d/g.c).groupby(level=cell).sum(); exp_ = (pi*a0.reindex(g.index)*g.m/m0.reindex(g.index)).groupby(level=cell).sum()
        rows.append(dict(year=yr, observed=1000*(w*obs).sum(), expected=1000*(w*exp_).sum()))
    return pd.DataFrame(rows).set_index('year')
chg = lambda s, a, b: (100*(s.loc[b, 'observed']/s.loc[a, 'observed'] - 1), 100*(s.loc[b, 'expected']/s.loc[a, 'observed'] - 1))
world = series(SR, 1990, an); can = series([101], 1990, an); hic = series([64], 1990, an)
# global deaths and cases equal the sums over super-regions
for col in ['d', 'c']:
    s_ = an[an.index.get_level_values('location_id').isin(SR)][col].groupby(level=['sex_id', 'age_id', 'year']).sum()
    g_ = an.xs(1, level='location_id')[col].reindex(s_.index); assert ((s_ - g_).abs()/g_).max() < 1e-6, col
fig = pd.concat({'World': world, 'High-income': hic, 'Canada': can}, names=['panel']).reset_index()
base_ = fig.groupby('panel')['observed'].transform('first'); fig['observed_idx'] = 100*fig['observed']/base_
fig['expected_idx'] = 100*fig['expected']/base_
# ---------------- 204 countries, 2000 and 2023 (2010 also)
R00, A00, R10 = rz(f'{gbd_annual}/IHME-GBD_2023_DATA-80a09abe-1.zip'), rz(f'{gbd_annual}/IHME-GBD_2023_DATA-532d5fea-1.zip'), rz(f'{gbd_annual}/IHME-GBD_2023_DATA-4af37d63-1.zip')
R23, A1923 = rz(f'{gbd_2019_2023}/IHME-GBD_2023_DATA-aaf6002d-1.zip'), rz(f'{gbd_2019_2023}/IHME-GBD_2023_DATA-8e54c028-1.zip')
lm = pd.read_csv(locmap, encoding='utf-8-sig'); countries = set(lm.location_id); srmap = lm.set_index('location_id').super_region_id
def cells(r, acx, yr):
    x = r[(r.metric_id == 1) & r.sex_id.isin([1, 2]) & r.age_id.isin(AGE40) & (r.year == yr) & r.location_id.isin(countries)].pivot_table(index=['location_id'] + cell, columns='measure_id', values='val').rename(columns={1: 'd', 5: 'c'})
    m = acx[acx.sex_id.isin([1, 2]) & acx.age_id.isin(AGE40) & (acx.year == yr) & (acx.metric_id == 3) & acx.location_id.isin(countries)].set_index(['location_id'] + cell).val.rename('m')/1e5
    return x[['d', 'c']].join(m, how='inner').assign(year=yr).set_index('year', append=True)
cty = pd.concat([cells(R00, A00, 2000), cells(R10, A00, 2010), cells(R23, A1923, 2023)])
assert cty.groupby(level=['location_id', 'year']).size().eq(24).all() and cty.index.get_level_values('location_id').nunique() == 204
rows = []
for loc in [1] + SR:
    units = sorted(countries) if loc == 1 else sorted(srmap[srmap == loc].index)
    s = series(units, 2000, cty); o, e = chg(s, 2000, 2023)
    rows.append(dict(location_id=loc, obs_2000_2023=o, exp_2000_2023=e))
cs = pd.DataFrame(rows).set_index('location_id')
per = []
for loc in sorted(countries):
    s = series([loc], 2000, cty); o, e = chg(s, 2000, 2023)
    per.append(dict(location_id=loc, location=lm.set_index('location_id').loc[loc, 'location_name'], super_region_id=srmap[loc], obs=o, exp=e,
                    obs_over_exp=100*((1 + o/100)/(1 + e/100) - 1), cases_2023=cty.xs((loc, 2023), level=['location_id', 'year']).c.sum()))
per = pd.DataFrame(per); per.to_csv(out/'G11_country_observed_expected_2000_2023.csv', index=False)
# ---------------- table S9
tab = []
for loc in [1] + SR + HI11:
    s = world if loc == 1 else series([loc], 1990, an)
    o, e = chg(s, 1990, 2023)
    if loc in cs.index: o2, e2 = cs.loc[loc, 'obs_2000_2023'], cs.loc[loc, 'exp_2000_2023']
    else: r_ = per.set_index('location_id').loc[loc]; o2, e2 = r_.obs, r_.exp
    tab.append(dict(location_id=loc, location=NAMES[loc], rate_1990=s.loc[1990, 'observed'], rate_2023=s.loc[2023, 'observed'], obs_1990_2023=o, exp_1990_2023=e,
                    obs_over_exp_1990_2023=100*((1 + o/100)/(1 + e/100) - 1), obs_2000_2023=o2, exp_2000_2023=e2, obs_over_exp_2000_2023=100*((1 + o2/100)/(1 + e2/100) - 1)))
tab = pd.DataFrame(tab); tab.to_csv(out/'G11_observed_expected_by_location.csv', index=False)
fd.mkdir(parents=True, exist_ok=True); fig.to_csv(fd/'Figure3_gbd_observed_expected.csv', index=False)
q = lambda x: [float(x.median()), float(x.quantile(.25)), float(x.quantile(.75))]
K = dict(world_1990_2023=chg(world, 1990, 2023), world_1990_2019=chg(world, 1990, 2019), canada_1990_2023=chg(can, 1990, 2023),
         world_countries_2000_2023=[float(cs.loc[1, 'obs_2000_2023']), float(cs.loc[1, 'exp_2000_2023'])],
         n_countries=int(len(per)), n_obs_below_exp=int((per.obs < per.exp).sum()), obs_over_exp_median_iqr=q(per.obs_over_exp),
         hi11_obs_over_exp_range=[float(tab[tab.location_id.isin(HI11)].obs_over_exp_1990_2023.min()), float(tab[tab.location_id.isin(HI11)].obs_over_exp_1990_2023.max())],
         world_deaths_1990=float(an.xs((1, 1990), level=['location_id', 'year']).d.sum()), world_deaths_2023=float(an.xs((1, 2023), level=['location_id', 'year']).d.sum()),
         world_cases_1990=float(an.xs((1, 1990), level=['location_id', 'year']).c.sum()), world_cases_2023=float(an.xs((1, 2023), level=['location_id', 'year']).c.sum()))
allc = an.xs(1, level='location_id').m.unstack(['sex_id', 'age_id']).mul(w, axis=1).sum(axis=1)*1000
K['world_allcause_change_2019_2021'] = float(100*(allc.loc[2021]/allc.loc[2019] - 1)); K['world_allcause_change_1990_2023'] = float(100*(allc.loc[2023]/allc.loc[1990] - 1))
K['world_observed_change_2019_2021'] = float(100*(world.loc[2021, 'observed']/world.loc[2019, 'observed'] - 1))
K['world_expected_change_2019_2021'] = float(100*(world.loc[2021, 'expected']/world.loc[2019, 'expected'] - 1))
K['high_income_1990_2023'] = chg(hic, 1990, 2023)
hs = an[an.index.get_level_values('location_id').isin(SR)].c.groupby(level=['year', 'location_id']).sum().unstack()
K['high_income_case_share_1990_2023'] = [float(hs.loc[1990, 64]/hs.loc[1990].sum()), float(hs.loc[2023, 64]/hs.loc[2023].sum())]
K['world_share_expected_log'] = float(np.log(1 + K['world_1990_2023'][1]/100)/np.log(1 + K['world_1990_2023'][0]/100))
json.dump(K, open(out/'G10_numbers.json', 'w'), indent=1, default=float)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20)
print(json.dumps(K, indent=1, default=lambda x: round(float(x), 2)))
print(tab.round(2).to_string(index=False))
print(per.groupby('super_region_id').apply(lambda g: pd.Series(dict(n=len(g), below=(g.obs < g.exp).sum(), med=g.obs_over_exp.median()))).round(1))
print(fig[fig.year.isin([1990, 2000, 2010, 2019, 2020, 2021, 2023])].round(2).to_string(index=False))
