"""GBD 2023 trend analyses. Args: gbd_annual_dir gbd_2019_2023_dir loc_map_csv who_csv std_csv t3_csv out_dir
G1: GBD RA deaths (ages >=20) vs WHO-registered M05-M06 deaths, annual, 11 high-income countries.
G2: GBD RA deaths at ages >=40, Canada, annual 1990-2023.
G3: GBD RA deaths per expected death among people with RA (cases x all-cause mortality), ages >=40, 2000/2010/2019/2023,
    global, super-regions, selected countries; age-sex-standardised decomposition of the change in RA deaths per case.
G5: GBD age-standardised RA death and prevalence rates, 1990-2023."""
import sys, os, io, json, zipfile
from pathlib import Path
import numpy as np, pandas as pd
gbd_annual, gbd_2019_2023, locmap, whocsv, stdcsv, t3csv, out = sys.argv[1:8]; out = Path(out); out.mkdir(parents=True, exist_ok=True)
def rz(path):
    z = zipfile.ZipFile(path); return pd.concat([pd.read_csv(io.BytesIO(z.read(n)), low_memory=False) for n in z.namelist() if n.endswith('.csv')])
F = {k: rz(f'{gbd_annual}/IHME-GBD_2023_DATA-{v}-1.zip') for k, v in dict(P1A='04c719ed', P1B='d1117f20', P1C='f5ad253b', P1D='0e1e6620', P2E='6a682da3', P2F='80a09abe', P2G='4af37d63', P2H='532d5fea').items()}
V = {k: rz(f'{gbd_2019_2023}/IHME-GBD_2023_DATA-{v}-1.zip') for k, v in dict(R23='aaf6002d', R19='c82898aa', A1923='8e54c028').items()}
AGE40 = [13, 14, 15, 16, 17, 18, 19, 20, 30, 31, 32, 235]
K = {}
# ---------- consistency checks
ra_annual = pd.concat([F['P1A'], F['P1B'], F['P1C']]).drop_duplicates(['measure_id', 'metric_id', 'location_id', 'sex_id', 'age_id', 'year'])
chk = ra_annual[(ra_annual.measure_id == 1) & (ra_annual.metric_id == 1) & (ra_annual.location_id == 1) & (ra_annual.sex_id == 3) & ra_annual.age_id.isin(AGE40) & ra_annual.year.isin([2000, 2010])].groupby('year').val.sum()
chk2 = pd.concat([F['P2F'], F['P2G']]); chk2 = chk2[(chk2.measure_id == 1) & (chk2.metric_id == 1) & (chk2.location_id == 1) & (chk2.sex_id == 3) & chk2.age_id.isin(AGE40)].groupby('year').val.sum()
K['check_global_ra_deaths_2000_2010'] = dict(annual_files=chk.to_dict(), year_files=chk2.to_dict())
# ---------- G1
who = pd.read_csv(whocsv)
c11 = {71: 'Australia', 101: 'Canada', 78: 'Denmark', 79: 'Finland', 80: 'France', 86: 'Italy', 67: 'Japan', 89: 'Netherlands', 93: 'Sweden', 95: 'United Kingdom', 102: 'United States of America'}
g20 = ra_annual[(ra_annual.measure_id == 1) & (ra_annual.metric_id == 1) & (ra_annual.age_id == 37) & (ra_annual.sex_id == 3) & ra_annual.location_id.isin(c11)]
g20 = g20.assign(country=g20.location_id.map(c11))[['country', 'year', 'val', 'lower', 'upper']]
m1 = g20.merge(who.rename(columns={'Year': 'year'})[['country', 'year', 'M05_M06']], on=['country', 'year'])
m1['ratio'] = m1.val/m1.M05_M06; m1.to_csv(out/'G1_gbd_vs_registered_annual.csv', index=False)
rows = []
for c, d in m1.groupby('country'):
    d = d[d.year.between(2000, 2019)].sort_values('year')
    X = np.c_[np.ones(len(d)), d.year - 2000]
    bg = np.linalg.lstsq(X, np.log(d.val), rcond=None)[0][1]; br = np.linalg.lstsq(X, np.log(d.M05_M06), rcond=None)[0][1]
    dg, dr = np.diff(np.log(d.val.values)), np.diff(np.log(d.M05_M06.values))
    rows.append(dict(country=c, years=f'{d.year.min()}-{d.year.max()}', n=len(d), ratio_median=d.ratio.median(), ratio_min=d.ratio.min(), ratio_max=d.ratio.max(),
                     apc_gbd=100*(np.exp(bg) - 1), apc_registered=100*(np.exp(br) - 1), apc_diff=100*(np.exp(bg) - np.exp(br)),
                     r_levels=np.corrcoef(d.val, d.M05_M06)[0, 1], r_changes=np.corrcoef(dg, dr)[0, 1] if len(dg) > 2 else np.nan))
g1 = pd.DataFrame(rows); g1.to_csv(out/'G1_summary.csv', index=False)
# ---------- G2
can = ra_annual[(ra_annual.measure_id == 1) & (ra_annual.metric_id == 1) & (ra_annual.location_id == 101) & (ra_annual.sex_id == 3) & ra_annual.age_id.isin(AGE40)]
g2 = can.groupby('year')[['val', 'lower', 'upper']].sum().reset_index(); g2.to_csv(out/'G2_canada_gbd_ge40_annual.csv', index=False)
# ---------- G3
lm = pd.read_csv(locmap, encoding='utf-8-sig')[['location_id', 'super_region']]
def year_cells(yr):
    if yr in (2000, 2010):
        ra = pd.concat([F['P2F'], F['P2G']]); ra = ra[ra.year == yr]; ac = F['P2H']; ac = ac[ac.year == yr]
    else:
        ra = V['R19'] if yr == 2019 else V['R23']; ac = V['A1923']; ac = ac[ac.year == yr]
    ra = ra[(ra.metric_id == 1) & ra.sex_id.isin([1, 2]) & ra.age_id.isin(AGE40)]
    p = ra.pivot_table(index=['location_id', 'sex_id', 'age_id'], columns='measure_id', values='val').rename(columns={1: 'ra_deaths', 5: 'cases', 3: 'yld', 4: 'yll'})
    a = ac[ac.sex_id.isin([1, 2]) & ac.age_id.isin(AGE40)].pivot_table(index=['location_id', 'sex_id', 'age_id'], columns='metric_id', values='val').rename(columns={1: 'all_deaths', 3: 'all_rate'})
    a['pop'] = a.all_deaths/(a.all_rate/1e5)
    c = p.join(a[['all_deaths', 'pop']], how='inner').reset_index(); c['year'] = yr
    return c
cells = pd.concat([year_cells(y) for y in [2000, 2010, 2019, 2023]])
assert cells.groupby('year').size().nunique() == 1, cells.groupby('year').size()
w = cells[(cells.year == 2019) & (cells.location_id == 1)].set_index(['sex_id', 'age_id']).cases; w = w/w.sum()
countries = cells[cells.location_id != 1].merge(lm, on='location_id', how='left'); assert countries.super_region.notna().all()
groups = {'Global': cells[cells.location_id == 1].assign(g='Global')}
sr = countries.groupby(['super_region', 'year', 'sex_id', 'age_id'])[['ra_deaths', 'cases', 'all_deaths', 'pop', 'yll', 'yld']].sum().reset_index().rename(columns={'super_region': 'g'})
sel = {101: 'Canada', 80: 'France', 78: 'Denmark', 86: 'Italy'}
cs = countries[countries.location_id.isin(sel)].assign(g=lambda d: d.location_id.map(sel))
hi11 = countries[countries.location_id.isin(c11)].groupby(['year', 'sex_id', 'age_id'])[['ra_deaths', 'cases', 'all_deaths', 'pop', 'yll', 'yld']].sum().reset_index().assign(g='11 high-income countries')
allg = pd.concat([groups['Global'], sr, cs, hi11])
rows = []
for (gname, yr), d in allg.groupby(['g', 'year']):
    d = d.set_index(['sex_id', 'age_id']); d['m'] = d.all_deaths/d['pop']; E = (d.cases*d.m).sum()
    ww = w.reindex(d.index)
    std_att = 1000*(ww*d.ra_deaths/d.cases).sum(); std_bg = 1000*(ww*d.m).sum()
    rows.append(dict(group=gname, year=yr, cases=d.cases.sum(), ra_deaths=d.ra_deaths.sum(), expected=E, A_pct=100*d.ra_deaths.sum()/E,
                     crude_att_per_1000=1000*d.ra_deaths.sum()/d.cases.sum(), std_att_per_1000=std_att, std_bg_per_1000=std_bg, A_std_pct=100*std_att/std_bg,
                     share_ge80_cases=d[d.index.get_level_values('age_id').isin([30, 31, 32, 235])].cases.sum()/d.cases.sum()))
g3 = pd.DataFrame(rows)
def ch(gname, col, a=2000, b=2023):
    s = g3[g3.group == gname].set_index('year')[col]; return 100*(s[b]/s[a] - 1)
dec = []
for gname in g3.group.unique():
    att, bg, A = ch(gname, 'std_att_per_1000'), ch(gname, 'std_bg_per_1000'), ch(gname, 'A_std_pct')
    la, lb = np.log(1 + att/100), np.log(1 + bg/100)
    dec.append(dict(group=gname, pct_change_std_att=att, pct_change_std_bg=bg, pct_change_A_std=A, share_of_decline_from_background=lb/la if la < 0 else np.nan,
                    pct_change_A_crude=ch(gname, 'A_pct'), pct_change_A_std_2010_23=ch(gname, 'A_std_pct', 2010, 2023), pct_change_ra_deaths=ch(gname, 'ra_deaths'), pct_change_cases=ch(gname, 'cases')))
dec = pd.DataFrame(dec); g3.to_csv(out/'G3_attribution_ratio_by_year.csv', index=False); dec.to_csv(out/'G3_decomposition_2000_2023.csv', index=False)
# Canadian surveillance comparators
s = pd.read_csv(stdcsv).set_index('year'); t3 = pd.read_csv(t3csv); t3a = t3[t3.ages == 'All ≥40'].set_index('year')
obs = pd.DataFrame({'std_certified_per_1000': s.std_certified_per_1000, 'std_pop_mort_per_1000': t3a.population_rate, 'srr': t3a.srr})
obs['A_obs_std_pct'] = 100*obs.std_certified_per_1000/obs.std_pop_mort_per_1000; obs['srr_minus_1'] = obs.srr - 1
obs.to_csv(out/'G3_canada_surveillance_comparators.csv')
K['canada_obs_change_2010_23'] = dict(A_obs_std=100*(obs.A_obs_std_pct[2023]/obs.A_obs_std_pct[2010] - 1), srr_minus_1=100*(obs.srr_minus_1[2023]/obs.srr_minus_1[2010] - 1),
                                      std_certified=100*(obs.std_certified_per_1000[2023]/obs.std_certified_per_1000[2010] - 1), std_pop=100*(obs.std_pop_mort_per_1000[2023]/obs.std_pop_mort_per_1000[2010] - 1))
K['canada_obs_change_2007_23'] = dict(A_obs_std=100*(obs.A_obs_std_pct[2023]/obs.A_obs_std_pct[2007] - 1), srr_minus_1=100*(obs.srr_minus_1[2023]/obs.srr_minus_1[2007] - 1))
# ---------- G5
asr = F['P2E']; asr = asr[(asr.location_id == 1) & (asr.sex_id == 3)].pivot_table(index='year', columns='measure_id', values='val')
K['global_asr'] = {'deaths_1990': asr.loc[1990, 1], 'deaths_2023': asr.loc[2023, 1], 'pct_deaths': 100*(asr.loc[2023, 1]/asr.loc[1990, 1] - 1),
                   'prev_1990': asr.loc[1990, 5], 'prev_2023': asr.loc[2023, 5], 'pct_prev': 100*(asr.loc[2023, 5]/asr.loc[1990, 5] - 1),
                   'pct_deaths_2000_23': 100*(asr.loc[2023, 1]/asr.loc[2000, 1] - 1)}
glob20 = ra_annual[(ra_annual.measure_id == 1) & (ra_annual.metric_id == 1) & (ra_annual.location_id == 1) & (ra_annual.sex_id == 3) & (ra_annual.age_id == 37)].set_index('year')
K['global_deaths_20plus'] = {int(y): [glob20.loc[y, 'val'], glob20.loc[y, 'lower'], glob20.loc[y, 'upper']] for y in [1990, 2000, 2010, 2019, 2023]}
json.dump(K, open(out/'gbd_trend_numbers.json', 'w'), indent=1, default=float)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(json.dumps(K, indent=1, default=float))
print(g1.round(3).to_string(index=False))
print(g3.round(3).to_string(index=False))
print(dec.round(2).to_string(index=False))
print(obs.round(3).loc[[2007, 2010, 2019, 2023]].to_string())
print(g2[g2.year >= 2007].round(1).to_string(index=False))
