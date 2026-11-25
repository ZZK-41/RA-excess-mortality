"""Annual 1990-2023 GBD RA deaths per person with RA vs all-cause mortality (ages >=40, fixed weights), and Canada GBD per-case standardised
to CCDSS 2007-09 weights. Args: gbd_annual_dir gbd_2019_2023_dir ccdss_raw out_dir"""
import sys, io, json, zipfile
from pathlib import Path
import numpy as np, pandas as pd
gbd_annual, gbd_2019_2023, ccdss, out = sys.argv[1:5]; out = Path(out); out.mkdir(parents=True, exist_ok=True)
def rz(p):
    z = zipfile.ZipFile(p); return pd.concat([pd.read_csv(io.BytesIO(z.read(n)), low_memory=False) for n in z.namelist() if n.endswith('.csv')])
ra = pd.concat([rz(f'{gbd_annual}/IHME-GBD_2023_DATA-{x}-1.zip') for x in ['04c719ed', 'd1117f20', 'f5ad253b']]).drop_duplicates(['measure_id', 'metric_id', 'location_id', 'sex_id', 'age_id', 'year'])
pv = pd.concat([rz(f'{gbd_annual}/IHME-GBD_2023_DATA-d9e5f203-1.zip'), rz(f'{gbd_annual}/IHME-GBD_2023_DATA-a4d3c489-1.zip')]); pv = pv[pv.measure_id == 5]; ac = rz(f'{gbd_annual}/IHME-GBD_2023_DATA-0e1e6620-1.zip')
P2F, P2G, R19, R23 = [rz(p) for p in [f'{gbd_annual}/IHME-GBD_2023_DATA-80a09abe-1.zip', f'{gbd_annual}/IHME-GBD_2023_DATA-4af37d63-1.zip', f'{gbd_2019_2023}/IHME-GBD_2023_DATA-c82898aa-1.zip', f'{gbd_2019_2023}/IHME-GBD_2023_DATA-aaf6002d-1.zip']]
AGE40 = [13, 14, 15, 16, 17, 18, 19, 20, 30, 31, 32, 235]; SR = [4, 31, 64, 103, 137, 158, 166]
NAMES = {1: 'Global', 64: 'High-income', 31: 'Central Europe, Eastern Europe, and Central Asia', 103: 'Latin America and Caribbean', 137: 'North Africa and Middle East',
         158: 'South Asia', 4: 'Southeast Asia, East Asia, and Oceania', 166: 'Sub-Saharan Africa', 101: 'Canada', 78: 'Denmark', 80: 'France', 86: 'Italy',
         71: 'Australia', 79: 'Finland', 67: 'Japan', 89: 'Netherlands', 93: 'Sweden', 95: 'United Kingdom', 102: 'United States of America'}
key = ['location_id', 'sex_id', 'age_id', 'year']
D = ra[(ra.measure_id == 1) & (ra.metric_id == 1) & ra.sex_id.isin([1, 2]) & ra.age_id.isin(AGE40)].set_index(key).val.rename('d')
C = pv[(pv.metric_id == 1) & pv.sex_id.isin([1, 2]) & pv.age_id.isin(AGE40)].set_index(key).val.rename('c')
Cg = C[C.index.get_level_values('location_id').isin(SR)].groupby(level=['sex_id', 'age_id', 'year']).sum()
Cglob = C.xs(1, level='location_id'); assert ((Cg - Cglob.reindex(Cg.index)).abs()/Cglob.reindex(Cg.index)).max() < 1e-9, 'global prevalence differs from sum of super-regions'
chk = {}
for yr, src in [(2000, P2F), (2010, P2G), (2019, R19), (2023, R23)]:
    s = src[(src.measure_id == 5) & (src.metric_id == 1) & (src.location_id == 1) & src.sex_id.isin([1, 2]) & src.age_id.isin(AGE40)].val.sum()
    chk[yr] = dict(sum_of_super_regions=float(Cg.xs(yr, level='year').sum()), global_annual_file=float(Cglob.xs(yr, level='year').sum()), global_year_file=float(s))
print('global prevalence check', chk)
A = ac[ac.sex_id.isin([1, 2]) & ac.age_id.isin(AGE40)].pivot_table(index=key, columns='metric_id', values='val'); M = (A[3]/1e5).rename('m')
df = pd.concat([D, C, M], axis=1).dropna().reset_index()
assert df.groupby(['location_id', 'year']).size().eq(24).all(), df.groupby(['location_id', 'year']).size().value_counts()
w = df[(df.location_id == 1) & (df.year == 2019)].set_index(['sex_id', 'age_id']).c; w = w/w.sum()
rows = []
for (loc, yr), g in df.groupby(['location_id', 'year']):
    g = g.set_index(['sex_id', 'age_id']); ww = w.reindex(g.index)
    att = 1000*(ww*g.d/g.c).sum(); bg = 1000*(ww*g.m).sum()
    rows.append(dict(location_id=loc, location=NAMES[loc], year=yr, ra_deaths=g.d.sum(), cases=g.c.sum(), std_att=att, std_bg=bg, ratio=att/bg))
s = pd.DataFrame(rows)
for col in ['std_att', 'std_bg', 'ratio']: s[f'{col}_idx'] = 100*s[col]/s.groupby('location_id')[col].transform(lambda x: x.iloc[0])
s.to_csv(out/'G6_annual_std_1990_2023.csv', index=False)
ch = []
for loc, g in s.groupby('location_id'):
    g = g.set_index('year'); r = dict(location=NAMES[loc])
    for a, b in [(1990, 2023), (2000, 2023), (2010, 2023), (1990, 2010)]:
        for col in ['std_att', 'std_bg', 'ratio']: r[f'{col}_{a}_{b}'] = 100*(g.loc[b, col]/g.loc[a, col] - 1)
    ch.append(r)
ch = pd.DataFrame(ch); ch.to_csv(out/'G6_changes.csv', index=False)
# Canada: GBD per-case RA deaths standardised to CCDSS 2007-09 RA weights (ages 40-44 ... 85-89, 90+)
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda x: pd.to_numeric(x.astype(str).str.replace(',', ''), errors='coerce')
t['n1'] = nm(t['Population of mortality with the condition'])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[c5.age >= 40]
c5['sex_id'] = c5.Sex.map({'Males': 1, 'Females': 2})
wc = c5[c5.year.between(2007, 2009)].groupby(['sex_id', 'age']).n1.sum(); wc = wc/wc.sum()
AMAP = {13: 40, 14: 45, 15: 50, 16: 55, 17: 60, 18: 65, 19: 70, 20: 75, 30: 80, 31: 85, 32: 90, 235: 90}
cg = df[df.location_id == 101].assign(age=lambda d: d.age_id.map(AMAP)).groupby(['year', 'sex_id', 'age'])[['d', 'c']].sum().reset_index()
cs = []
for yr, g in cg.groupby('year'):
    g = g.set_index(['sex_id', 'age']); cs.append(dict(year=yr, gbd_std_per_1000=1000*(wc.reindex(g.index)*g.d/g.c).sum(), gbd_deaths=g.d.sum(), gbd_cases=g.c.sum()))
cs = pd.DataFrame(cs); base = cs[cs.year.between(2007, 2009)].gbd_std_per_1000.mean(); cs['gbd_index_0709'] = 100*cs.gbd_std_per_1000/base
cs.to_csv(out/'G7_canada_gbd_ccdss_weights.csv', index=False)
pm = lambda a, b: 100*(cs[cs.year.isin(b)].gbd_std_per_1000.mean()/cs[cs.year.isin(a)].gbd_std_per_1000.mean() - 1)
K = dict(global_prev_check=chk, canada_gbd_change_0709_2123=pm([2007, 2008, 2009], [2021, 2022, 2023]), canada_gbd_change_0709_1719=pm([2007, 2008, 2009], [2017, 2018, 2019]),
         canada_gbd_change_1719_2123=pm([2017, 2018, 2019], [2021, 2022, 2023]))
json.dump(K, open(out/'G6_numbers.json', 'w'), indent=1, default=float)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 20)
print(json.dumps({k: v for k, v in K.items() if k != 'global_prev_check'}, indent=1, default=float))
print(ch[['location', 'std_att_1990_2023', 'std_bg_1990_2023', 'ratio_1990_2023', 'std_att_2000_2023', 'std_bg_2000_2023', 'ratio_2000_2023', 'ratio_2010_2023', 'ratio_1990_2010']].round(1).to_string(index=False))
print(s[s.location_id == 1][['year', 'ra_deaths', 'cases', 'std_att', 'std_bg', 'ratio']].iloc[::4].round(3).to_string(index=False))
