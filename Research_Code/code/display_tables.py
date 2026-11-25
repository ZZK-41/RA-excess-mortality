"""Numbers for the text and tables (display_numbers.json) and CIs of the ratios for the Canadian table and the European comparisons.
Ratio intervals use a fixed seed for every call, so identical inputs give identical intervals.
Args: results_final_dir ccdss_raw statcan_zip gbd_dir out_dir"""
import sys, io, json, zipfile
from pathlib import Path
import numpy as np, pandas as pd
res, ccdss, statcan, gbddir, out = sys.argv[1:6]
res = Path(res); out = Path(out); out.mkdir(parents=True, exist_ok=True)
Z = 1.959964; B = 20000; SEED = 20261123
def ratio_ci(D, E, C):
    rng = np.random.default_rng(SEED); d = rng.poisson(D, B); c = np.maximum(rng.poisson(C, B), 1); r = (d - E)/c
    return float(np.quantile(r, .025)), float(np.quantile(r, .975))
def ratio_ci_general(ex, ex_lo, ex_hi, c, c_lo=None, c_hi=None, poisson_c=False):
    rng = np.random.default_rng(SEED); e = rng.normal(ex, (ex_hi - ex_lo)/(2*Z), B)
    if poisson_c: cc = np.maximum(rng.poisson(c, B), 1)
    else: cc = np.exp(rng.normal(np.log(c), (np.log(c_hi) - np.log(c_lo))/(2*Z), B))
    r = e/cc; return float(np.quantile(r, .025)), float(np.quantile(r, .975))
K = {}
# ---------------- Canada annual series (Figure 1, Table S1)
t1 = pd.read_csv(res/'T1_canada_annual_excess_vs_coded.csv')
ci = t1.apply(lambda r: ratio_ci(r.deaths, r.expected, r.coded_statcan), axis=1, result_type='expand'); t1['ratio_lo'], t1['ratio_hi'] = ci[0], ci[1]
t1['share_pct'] = 100*t1.coded_statcan/t1.excess; t1['share_lo'] = 100/t1.ratio_hi; t1['share_hi'] = 100/t1.ratio_lo
K['cumulative'] = dict(excess=float(t1.excess.sum()), coded=float(t1.coded_statcan.sum()), ratio=float(t1.excess.sum()/t1.coded_statcan.sum()),
                       share_pct=float(100*t1.coded_statcan.sum()/t1.excess.sum()))
K['annual_2023'] = t1[t1.year == 2023].iloc[0].to_dict(); K['annual_2019'] = t1[t1.year == 2019].iloc[0].to_dict()
K['ratio_range'] = [float(t1.ratio.min()), int(t1.loc[t1.ratio.idxmin(), 'year']), float(t1.ratio.max()), int(t1.loc[t1.ratio.idxmax(), 'year'])]
# trends in certified share: whole period, pre-pandemic, constant jurisdiction coverage (2010-23)
def share_trend(d):
    share = d.coded_statcan/d.excess; v = 1/d.coded_statcan + d.deaths/d.excess**2
    W = 1/v.values; X = np.c_[np.ones(len(d)), d.year.values - d.year.min()]; y = np.log(share.values)
    A = X.T@(X*W[:, None]); b = np.linalg.solve(A, X.T@(W*y)); r = y - X@b; disp = max(1, (W*r**2).sum()/(len(d) - 2)); se = np.sqrt(np.linalg.inv(A)[1, 1]*disp)
    return [100*(np.exp(b[1]) - 1), 100*(np.exp(b[1] - Z*se) - 1), 100*(np.exp(b[1] + Z*se) - 1)]
K['share_trend'] = {'2007-23': share_trend(t1), '2007-19': share_trend(t1[t1.year <= 2019]), '2010-23': share_trend(t1[t1.year >= 2010]), '2010-19': share_trend(t1[t1.year.between(2010, 2019)])}
for lab, yrs in [('2007-09', [2007, 2008, 2009]), ('2010-12', [2010, 2011, 2012]), ('2017-19', [2017, 2018, 2019]), ('2021-23', [2021, 2022, 2023])]:
    s = t1[t1.year.isin(yrs)]; K[f'share_{lab}'] = float(100*s.coded_statcan.sum()/s.excess.sum())
K['coded_per_1000'] = {int(y): float(1000*c/n) for y, c, n in zip(t1.year, t1.coded_statcan, t1.people_ra)}
K['excess_per_1000'] = {int(y): float(v) for y, v in zip(t1.year, t1.excess_per_1000)}
# ---------------- Canada 2023 by age and sex (Table 1, Figure 2)
t2 = pd.read_csv(res/'T2_canada_by_age_sex_2019_2023.csv')
ci = t2.apply(lambda r: ratio_ci(r.deaths, r.expected, r.coded_statcan) if r.coded_statcan > 0 else (np.nan, np.nan), axis=1, result_type='expand')
t2['ratio_lo'], t2['ratio_hi'] = ci[0], ci[1]
cons = json.load(open(res/'conservative_5yr.json'))
K['conservative'] = cons
# ---------------- International (Table 2, Figure 2)
t5 = pd.read_csv(res/'T5_international.csv')
rows = []
for _, r in t5.iterrows():
    if 'Veneto' in r.setting: lo, hi = ratio_ci_general(r.excess, r.excess_lo, r.excess_hi, r.coded, poisson_c=True)
    else: lo, hi = ratio_ci_general(r.excess, r.excess_lo, r.excess_hi, r.coded, r.coded_lo, r.coded_hi)
    rows.append(r.to_dict() | dict(ratio_lo_combined=lo, ratio_hi_combined=hi))
t5 = pd.DataFrame(rows)
# ---------------- Global (Figure 3, Table 2)
t6 = pd.read_csv(res/'T6_global_45_89.csv'); t6['gbd_yll_share'] = t6.gbd_yll/t6.daly_gbd
K['global_2023'] = t6[t6.year == 2023].set_index('scenario')[['cases', 'gbd_coded', 'modelled_excess', 'ratio', 'gbd_yll', 'gbd_yld', 'excess_yll', 'daly_gbd', 'daly_with_excess', 'daly_ratio', 'yll_share_with_excess', 'gbd_yll_share', 'share_80_89']].to_dict('index')
K['global_2019'] = t6[t6.year == 2019].set_index('scenario')[['modelled_excess', 'ratio', 'daly_ratio']].to_dict('index')
# ---------------- conservative Canadian YLL (excess per CCDSS case x GBD cases, by 5-year age and sex, x GBD YLL per death)
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
for a, b in [('d1', 'Counts of mortality with the condition'), ('n1', 'Population of mortality with the condition'), ('d0', 'Counts of mortality without the condition'), ('n0', 'Population of mortality without the condition')]: t[a] = nm(t[b])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[c5.age >= 40]
c5['sex'] = c5.Sex.map({'Females': 'Female', 'Males': 'Male'}); c5['ex'] = c5.d1 - c5.n1*(c5.d1 + c5.d0)/(c5.n1 + c5.n0)
z = zipfile.ZipFile(f'{gbddir}/IHME-GBD_2023_DATA-aaf6002d-1.zip'); m = [m for m in z.infolist() if m.filename.endswith('.csv')][0]
g = pd.read_csv(io.BytesIO(z.read(m)), low_memory=False)
g = g[(g.location_name == 'Canada') & (g.metric_name == 'Number') & g.sex_name.isin(['Female', 'Male'])]
FIVE = {f'{a}-{a+4} years': a for a in range(20, 95, 5)} | {'95+ years': 95}
g = g[g.age_name.isin(FIVE)].copy(); g['age'] = g.age_name.map(FIVE).clip(upper=90)
p = g.pivot_table(index=['sex_name', 'age'], columns='measure_name', values='val', aggfunc='sum').reset_index().rename(columns={'sex_name': 'sex'})
p['ypd'] = p['YLLs (Years of Life Lost)']/p['Deaths']
s = c5[c5.year == 2023].merge(p[['sex', 'age', 'ypd', 'Prevalence']], on=['sex', 'age'])
yll = json.load(open(res/'canada_yll_trend.json'))
K['canada_yll'] = yll | dict(conservative_excess_yll=float((s.ex/s.n1*s.Prevalence*s.ypd).sum()))
K['canada_yll']['conservative_daly_ratio'] = (yll['gbd_yld_20plus'] + K['canada_yll']['conservative_excess_yll'])/yll['gbd_daly_20plus']
K['canada_yll']['excess_yll_per_excess_death'] = yll['excess_yll_2023']/yll['excess_deaths_2023']
K['canada_yll']['gbd_yll_per_gbd_death'] = yll['gbd_yll_20plus']/yll['gbd_deaths_20plus']
# share of excess and deaths at >=80
s23 = c5[c5.year == 2023]; K['share_excess_ge80_2023'] = float(s23[s23.age >= 80].ex.sum()/s23.ex.sum()); K['share_people_ge80_2023'] = float(s23[s23.age >= 80].n1.sum()/s23.n1.sum())
K['share_excess_ge90_2023'] = float(s23[s23.age >= 90].ex.sum()/s23.ex.sum())
json.dump(K, open(out/'display_numbers.json', 'w'), indent=1, default=float)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(t1[['year', 'excess', 'coded_statcan', 'ratio', 'ratio_lo', 'ratio_hi', 'share_pct', 'share_lo', 'share_hi']].round(3).to_string(index=False))
print(t2[t2.year == 2023][['sex', 'ages', 'people_ra', 'deaths', 'expected', 'excess', 'excess_lo', 'excess_hi', 'smr', 'smr_lo', 'smr_hi', 'coded_statcan', 'ratio', 'ratio_lo', 'ratio_hi']].round(3).to_string(index=False))
print(t5[['setting', 'excess', 'coded', 'ratio', 'ratio_lo_combined', 'ratio_hi_combined']].round(3).to_string(index=False))
print(json.dumps({k: v for k, v in K.items() if k not in ('conservative', 'annual_2023', 'annual_2019')}, indent=1, default=float))
