"""95% CI for the conservative Canadian scenario (excess per person by 5-year age and sex x GBD case count), 2023 and 2019.
Args: ccdss_raw statcan_zip gbd_dir"""
import sys, io, zipfile
import numpy as np, pandas as pd
ccdss, statcan, gbddir = sys.argv[1:4]; B = 20000
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
for a, b in [('d1', 'Counts of mortality with the condition'), ('n1', 'Population of mortality with the condition'), ('d0', 'Counts of mortality without the condition'), ('n0', 'Population of mortality without the condition')]: t[a] = nm(t[b])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[c5.age >= 40]
c5['sex'] = c5.Sex.map({'Females': 'Female', 'Males': 'Male'}); c5['mpop'] = (c5.d1 + c5.d0)/(c5.n1 + c5.n0)
def gz(name):
    z = zipfile.ZipFile(f'{gbddir}/{name}'); m = [m for m in z.infolist() if m.filename.endswith('.csv')][0]
    return pd.read_csv(io.BytesIO(z.read(m)), low_memory=False)
g = pd.concat([gz('IHME-GBD_2023_DATA-aaf6002d-1.zip'), gz('IHME-GBD_2023_DATA-c82898aa-1.zip')])
g = g[(g.location_name == 'Canada') & (g.metric_name == 'Number') & g.sex_name.isin(['Female', 'Male']) & (g.measure_name == 'Prevalence')]
FIVE = {f'{a}-{a+4} years': a for a in range(20, 95, 5)} | {'95+ years': 95}
g = g[g.age_name.isin(FIVE)].copy(); g['age'] = g.age_name.map(FIVE).clip(upper=90)
sc = pd.read_csv(io.BytesIO(zipfile.ZipFile(statcan).read('13100150.csv')), low_memory=False)
sc = sc[sc['Cause of death (ICD-10)'].isin(['Seropositive rheumatoid arthritis [M05]', 'Other rheumatoid arthritis [M06]']) & sc.Sex.eq('Both sexes')]
amap = {f'{a} to {a+4} years': a for a in range(40, 90, 5)} | {'90 years and over': 90}
sc['age'] = sc['Age group'].map(amap); sc = sc[sc.age.notna()]
for yr in [2019, 2023]:
    p = g[g.year == yr].groupby(['sex_name', 'age']).val.sum().rename('G').reset_index().rename(columns={'sex_name': 'sex'})
    s = c5[c5.year == yr].merge(p, on=['sex', 'age']); C = float(sc[sc.REF_DATE == yr].VALUE.sum())
    rng = np.random.default_rng(20261123)
    D = rng.poisson(s.d1.values, (B, len(s))); f = s.G.values/s.n1.values
    cons = ((D - (s.n1*s.mpop).values)*f).sum(1); cc = np.maximum(rng.poisson(C, B), 1)
    point = float(((s.d1 - s.n1*s.mpop)*f).sum())
    r = cons/cc
    print(yr, 'conservative excess %.0f (%.0f-%.0f); certified %d; ratio %.2f (%.2f-%.2f); share %.1f%%' % (point, *np.quantile(cons, [.025, .975]), C, point/C, *np.quantile(r, [.025, .975]), 100*C/point))
