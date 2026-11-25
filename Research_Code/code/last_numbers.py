"""Remaining numbers: pooled 2007-23 ratio with 95% CI; conservative (GBD case count) ratios by sex for 2019 and 2023;
SRR 2010 vs 2019 standardised rates. Args: ccdss_raw statcan_zip gbd_dir t3_csv out_json"""
import sys, io, json, zipfile
import numpy as np, pandas as pd
ccdss, statcan, gbddir, t3, out = sys.argv[1:6]; B = 20000
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
for a, b in [('d1', 'Counts of mortality with the condition'), ('n1', 'Population of mortality with the condition'), ('d0', 'Counts of mortality without the condition'), ('n0', 'Population of mortality without the condition')]: t[a] = nm(t[b])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[c5.age >= 40]
c5['sex'] = c5.Sex.map({'Females': 'Female', 'Males': 'Male'}); c5['mpop'] = (c5.d1 + c5.d0)/(c5.n1 + c5.n0); c5['E'] = c5.n1*c5.mpop
sc = pd.read_csv(io.BytesIO(zipfile.ZipFile(statcan).read('13100150.csv')), low_memory=False)
sc = sc[sc['Cause of death (ICD-10)'].isin(['Seropositive rheumatoid arthritis [M05]', 'Other rheumatoid arthritis [M06]'])]
amap = {f'{a} to {a+4} years': a for a in range(40, 90, 5)} | {'90 years and over': 90}
sc['age'] = sc['Age group'].map(amap); sc['sex'] = sc.Sex.map({'Females': 'Female', 'Males': 'Male', 'Both sexes': 'Both'}); sc = sc[sc.age.notna()]
K = {}
rng = np.random.default_rng(20261123)
D, E, C = c5.d1.sum(), c5.E.sum(), sc[(sc.sex == 'Both') & sc.REF_DATE.between(2007, 2023)].VALUE.sum()
sim = (rng.poisson(D, B) - E)/np.maximum(rng.poisson(C, B), 1)
K['pooled_2007_23'] = dict(deaths=float(D), expected=float(E), excess=float(D - E), certified=float(C), ratio=float((D - E)/C), lo=float(np.quantile(sim, .025)), hi=float(np.quantile(sim, .975)), smr=float(D/E))
def gz(name):
    z = zipfile.ZipFile(f'{gbddir}/{name}'); m = [m for m in z.infolist() if m.filename.endswith('.csv')][0]
    return pd.read_csv(io.BytesIO(z.read(m)), low_memory=False)
g = pd.concat([gz('IHME-GBD_2023_DATA-aaf6002d-1.zip'), gz('IHME-GBD_2023_DATA-c82898aa-1.zip')])
g = g[(g.location_name == 'Canada') & (g.metric_name == 'Number') & g.sex_name.isin(['Female', 'Male']) & (g.measure_name == 'Prevalence')]
FIVE = {f'{a}-{a+4} years': a for a in range(20, 95, 5)} | {'95+ years': 95}
g = g[g.age_name.isin(FIVE)].copy(); g['age'] = g.age_name.map(FIVE).clip(upper=90)
for yr in [2019, 2023]:
    p = g[g.year == yr].groupby(['sex_name', 'age']).val.sum().rename('G').reset_index().rename(columns={'sex_name': 'sex'})
    s = c5[c5.year == yr].merge(p, on=['sex', 'age'])
    for sx in ['Female', 'Male']:
        q = s[s.sex == sx]; cons = float(((q.d1 - q.E)*q.G/q.n1).sum()); cc = float(sc[(sc.REF_DATE == yr) & (sc.sex == sx)].VALUE.sum())
        K[f'conservative_{yr}_{sx}'] = dict(conservative_excess=cons, certified=cc, ratio=cons/cc, full_excess=float((q.d1 - q.E).sum()), full_ratio=float((q.d1 - q.E).sum()/cc),
                                         ccdss_people=float(q.n1.sum()), gbd_people=float(q.G.sum()))
    q = s[s.age >= 80]; K[f'conservative_{yr}_share_ge80'] = float(((q.d1 - q.E)*q.G/q.n1).sum()/((s.d1 - s.E)*s.G/s.n1).sum())
    K[f'gbd_people_share_ge80_{yr}'] = float(s[s.age >= 80].G.sum()/s.G.sum())
s3 = pd.read_csv(t3); a = s3[s3.ages == 'All ≥40'].set_index('year')
K['std_rates'] = {int(y): dict(patient=float(a.loc[y, 'patient_rate']), population=float(a.loc[y, 'population_rate']), srr=float(a.loc[y, 'srr']), lo=float(a.loc[y, 'srr_lo']), hi=float(a.loc[y, 'srr_hi'])) for y in [2007, 2010, 2019, 2023]}
json.dump(K, open(out, 'w'), indent=1); print(json.dumps(K, indent=1))
