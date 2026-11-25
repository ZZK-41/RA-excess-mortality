"""Percentage change in standardised attributed and excess deaths per person with RA (Canada, ages >=40) between periods, with simulation CIs.
Args: ccdss_raw statcan_zip out_json"""
import sys, io, json, zipfile
import numpy as np, pandas as pd
ccdss, statcan, out = sys.argv[1:4]; rng = np.random.default_rng(20261123); B = 5000
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda x: pd.to_numeric(x.astype(str).str.replace(',', ''), errors='coerce')
for k, col in [('d1', 'Counts of mortality with the condition'), ('n1', 'Population of mortality with the condition'), ('d0', 'Counts of mortality without the condition'), ('n0', 'Population of mortality without the condition')]: t[k] = nm(t[col])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[c5.age >= 40].copy()
c5['sex'] = c5.Sex.map({'Females': 'F', 'Males': 'M'}); c5['N'] = c5.n1 + c5.n0
w = c5[c5.year.between(2007, 2009)].groupby(['age', 'sex']).n1.sum(); w = (w/w.sum()).rename('w').reset_index(); c5 = c5.merge(w, on=['age', 'sex'])
sc = pd.read_csv(io.BytesIO(zipfile.ZipFile(statcan).read('13100150.csv')), low_memory=False)
sc = sc[sc['Cause of death (ICD-10)'].isin(['Seropositive rheumatoid arthritis [M05]', 'Other rheumatoid arthritis [M06]']) & sc.Sex.isin(['Females', 'Males'])]
sc['age'] = sc['Age group'].map({f'{x} to {x+4} years': x for x in range(40, 90, 5)} | {'90 years and over': 90}); sc['sex'] = sc.Sex.map({'Females': 'F', 'Males': 'M'}); sc = sc[sc.age.notna()]
cc = sc.groupby(['REF_DATE', 'age', 'sex']).VALUE.sum().rename('C').reset_index().rename(columns={'REF_DATE': 'year'})
c5 = c5.merge(cc, on=['year', 'age', 'sex'], how='left').fillna({'C': 0})
def rates(g, D1, C, D0=None):
    D0 = g.d0.values if D0 is None else D0
    return 1000*(g.w*C/g.n1).sum(), 1000*(g.w*(D1/g.n1 - (D1 + D0)/g.N)).sum(), 1000*(g.w*(D1 + D0)/g.N).sum()
def period(yrs, sim=False):
    vals = []
    for _, g in c5[c5.year.isin(yrs)].groupby('year'):
        D1 = rng.poisson(g.d1.values) if sim else g.d1.values; C = rng.poisson(g.C.values) if sim else g.C.values
        D0 = rng.poisson(g.d0.values) if sim else g.d0.values
        vals.append(rates(g, D1, C, D0))
    return np.mean(vals, axis=0)
res = {}
for lab, a, b in [('2007-09_to_2017-19', [2007, 2008, 2009], [2017, 2018, 2019]), ('2007-09_to_2021-23', [2007, 2008, 2009], [2021, 2022, 2023]), ('2017-19_to_2021-23', [2017, 2018, 2019], [2021, 2022, 2023]), ('2007-09_to_2023', [2007, 2008, 2009], [2023])]:
    pa, pb = period(a), period(b); point = 100*(pb/pa - 1)
    sims = np.array([100*(period(b, True)/period(a, True) - 1) for _ in range(B)])
    lo, hi = np.quantile(sims, [.025, .975], axis=0)
    res[lab] = dict(certified=[point[0], lo[0], hi[0]], excess=[point[1], lo[1], hi[1]], whole_population=[point[2], lo[2], hi[2]])
    print(lab, 'certified %.1f%% (%.1f to %.1f); excess %.1f%% (%.1f to %.1f); whole-population mortality %.1f%% (%.1f to %.1f)' % (point[0], lo[0], hi[0], point[1], lo[1], hi[1], point[2], lo[2], hi[2]))
json.dump(res, open(out, 'w'), indent=1, default=float)
