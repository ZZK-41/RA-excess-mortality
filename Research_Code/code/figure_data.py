"""Canadian annual series and the ratios shown in figure 2 (excess deaths per attributed death, with 95% CIs). Deterministic (fixed seed per interval).
Args: ccdss_raw statcan_zip gbd_dir results_final_dir who_csv out_dir"""
import sys, io, json, zipfile
from pathlib import Path
import numpy as np, pandas as pd
ccdss, statcan, gbddir, resdir, whocsv, out = sys.argv[1:7]
resdir = Path(resdir); out = Path(out); out.mkdir(parents=True, exist_ok=True)
Z = 1.959964; B = 20000; SEED = 20261123
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
for a, b in [('d1', 'Counts of mortality with the condition'), ('n1', 'Population of mortality with the condition'), ('d0', 'Counts of mortality without the condition'), ('n0', 'Population of mortality without the condition')]: t[a] = nm(t[b])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[c5.age >= 40].copy()
c5['sex'] = c5.Sex.map({'Females': 'Female', 'Males': 'Male'}); c5['N'] = c5.n1 + c5.n0; c5['mpop'] = (c5.d1 + c5.d0)/c5.N; c5['E'] = c5.n1*c5.mpop
w = c5[c5.year.between(2007, 2009)].groupby(['age', 'sex']).n1.sum(); w = (w/w.sum()).rename('w').reset_index(); c5 = c5.merge(w, on=['age', 'sex'])
sc = pd.read_csv(io.BytesIO(zipfile.ZipFile(statcan).read('13100150.csv')), low_memory=False)
sc = sc[sc['Cause of death (ICD-10)'].isin(['Seropositive rheumatoid arthritis [M05]', 'Other rheumatoid arthritis [M06]'])]
amap = {f'{a} to {a+4} years': a for a in range(40, 90, 5)} | {'90 years and over': 90}
sc['age'] = sc['Age group'].map(amap); sc['sex'] = sc.Sex.map({'Females': 'Female', 'Males': 'Male', 'Both sexes': 'Both'}); sc = sc[sc.age.notna()]
cert = sc[sc.sex.isin(['Female', 'Male'])].groupby(['REF_DATE', 'age', 'sex']).VALUE.sum().rename('C').reset_index().rename(columns={'REF_DATE': 'year'})
c5 = c5.merge(cert, on=['year', 'age', 'sex'], how='left').fillna({'C': 0})
def ratio_ci(D, E, C):
    rng = np.random.default_rng(SEED); r = (rng.poisson(D, B) - E)/np.maximum(rng.poisson(C, B), 1)
    return float(np.quantile(r, .025)), float(np.quantile(r, .975))
def summ(d, label, group):
    D, E, C = d.d1.sum(), d.E.sum(), d.C.sum(); lo, hi = ratio_ci(D, E, C)
    return dict(group=group, label=label, deaths=D, expected=E, excess=D - E, excess_lo=D - E - Z*np.sqrt(D), excess_hi=D - E + Z*np.sqrt(D),
                smr=D/E, certified=C, ratio=(D - E)/C, lo=lo, hi=hi, share_pct=100*C/(D - E), share_lo=100/hi, share_hi=100/lo, kind='ci')
# ---- Figure 1: annual series
f1 = pd.DataFrame([summ(g, str(y), 'year') | dict(year=y, people_ra=g.n1.sum()) for y, g in c5.groupby('year')])
f1.to_csv(out/'Figure1_canada_annual.csv', index=False)
# ---- Figure S1: age-sex-standardised rates per 1000 people with RA (2007-09 weights), indexed to 2007-09
rows = []
for y, g in c5.groupby('year'):
    cr = 1000*(g.w*g.C/g.n1).sum(); vc = 1e6*(g.w**2*g.C/g.n1**2).sum()
    ex = 1000*(g.w*(g.d1/g.n1 - g.mpop)).sum(); ve = 1e6*(g.w**2*g.d1/g.n1**2).sum()
    rows.append(dict(year=y, certified_rate=cr, certified_lo=cr - Z*np.sqrt(vc), certified_hi=cr + Z*np.sqrt(vc), excess_rate=ex, excess_lo=ex - Z*np.sqrt(ve), excess_hi=ex + Z*np.sqrt(ve)))
s1 = pd.DataFrame(rows)
for k in ['certified', 'excess']:
    base = s1[s1.year.between(2007, 2009)][f'{k}_rate'].mean()
    for col in [f'{k}_rate', f'{k}_lo', f'{k}_hi']: s1[col.replace('rate', 'index') if col.endswith('rate') else col + '_index'] = 100*s1[col]/base
# ---- Figure 2: ratio of excess to certified deaths
f2 = [summ(c5, 'All ages ≥40 years', 'Canada, 2007–23 pooled')]
for lab, (lo_, hi_) in {'40–64 years': (40, 65), '65–79 years': (65, 80), '≥80 years': (80, 200)}.items():
    f2.append(summ(c5[(c5.age >= lo_) & (c5.age < hi_)], lab, 'Canada, 2007–23 pooled'))
f2 += [summ(c5[c5.sex == 'Female'], 'Women', 'Canada, 2007–23 pooled'), summ(c5[c5.sex == 'Male'], 'Men', 'Canada, 2007–23 pooled')]
f2.append(summ(c5[c5.year == 2023], 'All ages ≥40 years', 'Canada, 2023'))
# conservative: excess per person x GBD case counts by 5-year age and sex
def gz(name):
    z = zipfile.ZipFile(f'{gbddir}/{name}'); m = [m for m in z.infolist() if m.filename.endswith('.csv')][0]
    return pd.read_csv(io.BytesIO(z.read(m)), low_memory=False)
g = gz('IHME-GBD_2023_DATA-aaf6002d-1.zip')
g = g[(g.location_name == 'Canada') & (g.metric_name == 'Number') & g.sex_name.isin(['Female', 'Male']) & (g.measure_name == 'Prevalence') & (g.year == 2023)]
FIVE = {f'{a}-{a+4} years': a for a in range(20, 95, 5)} | {'95+ years': 95}
g = g[g.age_name.isin(FIVE)].copy(); g['age'] = g.age_name.map(FIVE).clip(upper=90)
p = g.groupby(['sex_name', 'age']).val.sum().rename('G').reset_index().rename(columns={'sex_name': 'sex'})
s = c5[c5.year == 2023].merge(p, on=['sex', 'age']); f = (s.G/s.n1).values; C23 = s.C.sum()
rng = np.random.default_rng(SEED); Dm = rng.poisson(s.d1.values, (B, len(s))); cons = ((Dm - s.E.values)*f).sum(1); r = cons/np.maximum(rng.poisson(C23, B), 1)
pt = float(((s.d1 - s.E)*f).sum())
f2.append(dict(group='Canada, 2023', label='Case numbers restricted to GBD estimate', excess=pt, excess_lo=float(np.quantile(cons, .025)), excess_hi=float(np.quantile(cons, .975)),
               certified=C23, ratio=pt/C23, lo=float(np.quantile(r, .025)), hi=float(np.quantile(r, .975)), share_pct=100*C23/pt, kind='ci'))
# Europe: registered underlying-cause deaths from the WHO Mortality Database
who = pd.read_csv(whocsv); t5 = pd.read_csv(resdir/'T5_international.csv')
fr_c = float(who[(who.country == 'France') & (who.Year == 2019)].M05_M06.iloc[0]); dk = who[(who.country == 'Denmark') & who.Year.between(2000, 2018)]
for _, rr in t5.iterrows():
    if 'France' in rr.setting: C, yrs, lab, grp = fr_c, 1, 'SNDS, 2019', 'France'
    elif 'M05' in rr.setting: C, yrs, lab, grp = float(dk.M05_M06.sum()), len(dk), 'Hospital-coded M05, 2000–18', 'Denmark'
    elif 'M06' in rr.setting: C, yrs, lab, grp = float(dk.M05_M06.sum()), len(dk), 'Hospital-coded M06, 2000–18', 'Denmark'
    else: C, yrs, lab, grp = float(rr.coded), 1, 'Prevalent cohort, 2012–20', 'Veneto, Italy'
    rng = np.random.default_rng(SEED); e = rng.normal(rr.excess, (rr.excess_hi - rr.excess_lo)/(2*Z), B); cc = np.maximum(rng.poisson(C, B), 1)/yrs
    q = e/cc; f2.append(dict(group=grp, label=lab, deaths=rr.deaths, expected=rr.expected, excess=rr.excess, excess_lo=rr.excess_lo, excess_hi=rr.excess_hi, smr=rr.smr,
                             certified=C/yrs, ratio=rr.excess/(C/yrs), lo=float(np.quantile(q, .025)), hi=float(np.quantile(q, .975)), share_pct=100*(C/yrs)/rr.excess,
                             gbd_deaths=rr.coded, ratio_vs_gbd=rr.ratio, kind='ci'))
t6 = pd.read_csv(resdir/'T6_global_45_89.csv'); g23 = t6[t6.year == 2023]
f2.append(dict(group='World, 2023 (modelled)', label='Ages 45–89 years, 204 locations', excess=g23.modelled_excess.median(), excess_lo=g23.modelled_excess.min(), excess_hi=g23.modelled_excess.max(),
               certified=g23.gbd_coded.iloc[0], ratio=np.nan, lo=g23.ratio.min(), hi=g23.ratio.max(), kind='range'))
f2 = pd.DataFrame(f2); f2.to_csv(out/'Figure2_ratios.csv', index=False)
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 30)
print(f2[['group', 'label', 'excess', 'certified', 'ratio', 'lo', 'hi', 'share_pct']].round(2).to_string(index=False))
print(s1.round(3).to_string(index=False))
