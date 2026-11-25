"""Age-sex-standardised certified RA death rate per 1000 people with RA vs standardised excess rate (fixed 2007-09 RA weights).
Args: ccdss_raw statcan_zip out_csv"""
import sys, io, zipfile
import numpy as np, pandas as pd
ccdss, statcan, out = sys.argv[1:4]; Z = 1.959964
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
for a, b in [('d1', 'Counts of mortality with the condition'), ('n1', 'Population of mortality with the condition'), ('d0', 'Counts of mortality without the condition'), ('n0', 'Population of mortality without the condition')]: t[a] = nm(t[b])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[c5.age >= 40].copy()
c5['sex'] = c5.Sex.map({'Females': 'Female', 'Males': 'Male'}); c5['N'] = c5.n1 + c5.n0; c5['mpop'] = (c5.d1 + c5.d0)/c5.N
w = c5[c5.year.between(2007, 2009)].groupby(['age', 'sex']).n1.sum(); w = (w/w.sum()).rename('w').reset_index(); c5 = c5.merge(w, on=['age', 'sex'])
sc = pd.read_csv(io.BytesIO(zipfile.ZipFile(statcan).read('13100150.csv')), low_memory=False)
sc = sc[sc['Cause of death (ICD-10)'].isin(['Seropositive rheumatoid arthritis [M05]', 'Other rheumatoid arthritis [M06]']) & sc.Sex.isin(['Females', 'Males'])]
amap = {f'{a} to {a+4} years': a for a in range(40, 90, 5)} | {'90 years and over': 90}
sc['age'] = sc['Age group'].map(amap); sc['sex'] = sc.Sex.map({'Females': 'Female', 'Males': 'Male'}); sc = sc[sc.age.notna()]
cc = sc.groupby(['REF_DATE', 'age', 'sex']).VALUE.sum().rename('C').reset_index().rename(columns={'REF_DATE': 'year'})
c5 = c5.merge(cc, on=['year', 'age', 'sex'], how='left').fillna({'C': 0})
rows = []
for yr, g in c5.groupby('year'):
    cert = (g.w*g.C/g.n1).sum()*1000; vc = (g.w**2*g.C/g.n1**2).sum()*1e6
    gap = (g.w*(g.d1/g.n1 - g.mpop)).sum()*1000; vg = (g.w**2*g.d1/g.n1**2).sum()*1e6
    rows.append(dict(year=yr, std_certified_per_1000=cert, cert_lo=cert - Z*np.sqrt(vc), cert_hi=cert + Z*np.sqrt(vc), std_excess_per_1000=gap,
                     gap_lo=gap - Z*np.sqrt(vg), gap_hi=gap + Z*np.sqrt(vg), std_share_pct=100*cert/gap))
r = pd.DataFrame(rows); r.to_csv(out, index=False); print(r.round(3).to_string(index=False))
def pct(a, b): return 100*(b/a - 1)
for col in ['std_certified_per_1000', 'std_excess_per_1000', 'std_share_pct']:
    a = r[r.year.between(2007, 2009)][col].mean(); b = r[r.year.between(2017, 2019)][col].mean(); c = r[r.year.between(2021, 2023)][col].mean()
    print(f'{col}: 2007-09 {a:.3f}  2017-19 {b:.3f} ({pct(a, b):+.1f}%)  2021-23 {c:.3f} ({pct(a, c):+.1f}%)')
for col in ['std_certified_per_1000', 'std_excess_per_1000']:
    d = r[r.year <= 2019]; X = np.c_[np.ones(len(d)), d.year - 2007]; y = np.log(d[col]); b = np.linalg.lstsq(X, y, rcond=None)[0]
    res = y - X@b; s2 = (res**2).sum()/(len(d) - 2); se = np.sqrt(s2*np.linalg.inv(X.T@X)[1, 1])
    print(f'{col} log-linear APC 2007-19: {100*(np.exp(b[1])-1):.2f}% ({100*(np.exp(b[1]-Z*se)-1):.2f} to {100*(np.exp(b[1]+Z*se)-1):.2f})')
