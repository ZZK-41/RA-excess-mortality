"""Canada: (1) years of life lost (YLL) for observed excess deaths, using GBD 2023 YLL per death by age and sex;
(2) trend in the share of excess deaths certified as RA, 2007-23.
Args: ccdss_raw statcan_zip gbd_dir out_json"""
import sys, io, json, zipfile
import numpy as np, pandas as pd
ccdss, statcan, gbddir, out = sys.argv[1:5]; Z = 1.959964
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
for a, b in [('d1', 'Counts of mortality with the condition'), ('n1', 'Population of mortality with the condition'), ('d0', 'Counts of mortality without the condition'), ('n0', 'Population of mortality without the condition')]: t[a] = nm(t[b])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[c5.age >= 40]
c5['sex'] = c5.Sex.map({'Females': 'Female', 'Males': 'Male'}); c5['ex'] = c5.d1 - c5.n1*(c5.d1 + c5.d0)/(c5.n1 + c5.n0)
def gz(name):
    z = zipfile.ZipFile(f'{gbddir}/{name}'); m = [m for m in z.infolist() if m.filename.endswith('.csv')][0]
    return pd.read_csv(io.BytesIO(z.read(m)), low_memory=False)
g = gz('IHME-GBD_2023_DATA-aaf6002d-1.zip')
g = g[(g.location_name == 'Canada') & (g.metric_name == 'Number') & g.sex_name.isin(['Female', 'Male'])]
FIVE = {f'{a}-{a+4} years': a for a in range(20, 95, 5)} | {'95+ years': 95}
g = g[g.age_name.isin(FIVE)].copy(); g['age'] = g.age_name.map(FIVE).clip(upper=90)
p = g.pivot_table(index=['sex_name', 'age'], columns='measure_name', values='val', aggfunc='sum').reset_index()
p['ypd'] = p['YLLs (Years of Life Lost)']/p['Deaths']
s = c5[c5.year == 2023].merge(p[['sex_name', 'age', 'ypd']].rename(columns={'sex_name': 'sex'}), on=['sex', 'age'])
res = dict(excess_deaths_2023=float(s.ex.sum()), excess_yll_2023=float((s.ex*s.ypd).sum()))
gg = gz('IHME-GBD_2023_DATA-aaf6002d-1.zip'); gg = gg[(gg.location_name == 'Canada') & (gg.metric_name == 'Number') & gg.sex_name.eq('Both') & gg.age_name.eq('20+ years')]
for msr, key in [('YLLs (Years of Life Lost)', 'gbd_yll_20plus'), ('YLDs (Years Lived with Disability)', 'gbd_yld_20plus'), ('Deaths', 'gbd_deaths_20plus')]:
    res[key] = float(gg[gg.measure_name == msr].val.iloc[0])
res['gbd_daly_20plus'] = res['gbd_yll_20plus'] + res['gbd_yld_20plus']
res['daly_with_excess'] = res['gbd_yld_20plus'] + res['excess_yll_2023']
res['daly_ratio'] = res['daly_with_excess']/res['gbd_daly_20plus']; res['yll_share_gbd'] = res['gbd_yll_20plus']/res['gbd_daly_20plus']
res['yll_share_with_excess'] = res['excess_yll_2023']/res['daly_with_excess']
# trend in certified share of excess deaths
sc = pd.read_csv(io.BytesIO(zipfile.ZipFile(statcan).read('13100150.csv')), low_memory=False)
sc = sc[sc['Cause of death (ICD-10)'].isin(['Seropositive rheumatoid arthritis [M05]', 'Other rheumatoid arthritis [M06]']) & sc.Sex.eq('Both sexes')]
amap = {f'{a} to {a+4} years': a for a in range(40, 90, 5)} | {'90 years and over': 90}
sc['age'] = sc['Age group'].map(amap); coded = sc[sc.age.notna()].groupby('REF_DATE').VALUE.sum()
a = c5.groupby('year').agg(ex=('ex', 'sum'), D=('d1', 'sum')); a['coded'] = coded.reindex(a.index)
a['share'] = a.coded/a.ex; a['var_log'] = 1/a.coded + a.D/a.ex**2
W = 1/a.var_log.values; X = np.c_[np.ones(len(a)), a.index.values - 2007]; y = np.log(a.share.values)
A = X.T@(X*W[:, None]); b = np.linalg.solve(A, X.T@(W*y)); r = y - X@b; disp = max(1, (W*r**2).sum()/(len(a) - 2)); se = np.sqrt(np.linalg.inv(A)[1, 1]*disp)
res['coded_share_apc_pct'] = [100*(np.exp(b[1]) - 1), 100*(np.exp(b[1] - Z*se) - 1), 100*(np.exp(b[1] + Z*se) - 1)]
res['coded_share_overdispersion'] = float(disp)
res['coded_per_1000_ra_2007'] = float(1000*a.coded.loc[2007]/c5[c5.year == 2007].n1.sum()); res['coded_per_1000_ra_2023'] = float(1000*a.coded.loc[2023]/c5[c5.year == 2023].n1.sum())
json.dump(res, open(out, 'w'), indent=1, default=float); print(json.dumps(res, indent=1, default=float))
