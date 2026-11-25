"""Conservative Canadian estimate: CCDSS excess deaths per person with RA, by 5-year age group and sex,
multiplied by GBD 2023 RA case counts in the same cells. Also GBD case counts for France and Denmark.
Args: ccdss_raw statcan_zip gbd_dir out_json"""
import sys, io, json, zipfile
import numpy as np, pandas as pd
ccdss, statcan, gbddir, out = sys.argv[1:5]
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
for a, b in [('d1', 'Counts of mortality with the condition'), ('n1', 'Population of mortality with the condition'), ('d0', 'Counts of mortality without the condition'), ('n0', 'Population of mortality without the condition')]: t[a] = nm(t[b])
can = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
print('CCDSS fiscal years with Canada 5-year cells:', sorted(can['Fiscal year'].unique())[:3], '...', sorted(can['Fiscal year'].unique())[-2:], len(can['Fiscal year'].unique()))
can['age'] = can['Age group'].str.extract(r'^(\d+)').astype(int); can['year'] = can['Fiscal year'].str[:4].astype(int)
print('earliest year with complete >=40 cells:', can[(can.age >= 40)].dropna(subset=['d1', 'n1', 'd0', 'n0']).groupby('year').size().loc[lambda s: s == 22].index.min())
c5 = can[can.age >= 40].copy(); c5['sex'] = c5.Sex.map({'Females': 'Female', 'Males': 'Male'})
c5['ex'] = c5.d1 - c5.n1*(c5.d1 + c5.d0)/(c5.n1 + c5.n0)
def gz(name):
    z = zipfile.ZipFile(f'{gbddir}/{name}'); m = [m for m in z.infolist() if m.filename.endswith('.csv')][0]
    return pd.read_csv(io.BytesIO(z.read(m)), low_memory=False)
g = pd.concat([gz('IHME-GBD_2023_DATA-aaf6002d-1.zip'), gz('IHME-GBD_2023_DATA-c82898aa-1.zip')]); g = g[g.metric_name == 'Number']
FIVE = {f'{a}-{a+4} years': a for a in range(20, 95, 5)} | {'95+ years': 95}
sc = pd.read_csv(io.BytesIO(zipfile.ZipFile(statcan).read('13100150.csv')), low_memory=False)
sc = sc[sc['Cause of death (ICD-10)'].isin(['Seropositive rheumatoid arthritis [M05]', 'Other rheumatoid arthritis [M06]']) & sc.Sex.eq('Both sexes')]
amap = {f'{a} to {a+4} years': a for a in range(0, 90, 5)} | {'90 years and over': 90}
sc['age'] = sc['Age group'].map(amap); sc = sc[sc.age.notna()]
res = {}
for yr in [2019, 2023]:
    gc = g[(g.location_name == 'Canada') & (g.year == yr) & g.sex_name.isin(['Female', 'Male']) & g.age_name.isin(FIVE)].copy()
    gc['age'] = gc.age_name.map(FIVE).clip(upper=90)
    p = gc.pivot_table(index=['sex_name', 'age'], columns='measure_name', values='val', aggfunc='sum').reset_index().rename(columns={'sex_name': 'sex'})
    s = c5[c5.year == yr].merge(p[['sex', 'age', 'Prevalence', 'Deaths']], on=['sex', 'age'], how='left')
    assert s.Prevalence.notna().all()
    s['cons'] = s.ex/s.n1*s.Prevalence
    for lab, lo, hi in [('ge40', 40, 200), ('45_89', 45, 90), ('40_64', 40, 65), ('65_79', 65, 80), ('ge80', 80, 200)]:
        q = s[(s.age >= lo) & (s.age < hi)]; coded = float(sc[(sc.REF_DATE == yr) & (sc.age >= lo) & (sc.age < hi)].VALUE.sum())
        res[f'{yr}_{lab}'] = dict(excess_ccdss=float(q.ex.sum()), people_ccdss=float(q.n1.sum()), people_gbd=float(q.Prevalence.sum()), gbd_deaths=float(q.Deaths.sum()),
                                  conservative_excess=float(q.cons.sum()), coded=coded, ratio_full=float(q.ex.sum()/coded), ratio_conservative=float(q.cons.sum()/coded),
                                  ratio_conservative_vs_gbd=float(q.cons.sum()/q.Deaths.sum()), case_ratio_ccdss_to_gbd=float(q.n1.sum()/q.Prevalence.sum()))
# GBD case counts and deaths for France and Denmark (2019, ages 20+), global YLL share at 45-89 is in T6
for loc in ['France', 'Denmark', 'Canada']:
    for yr in [2019, 2023]:
        d = g[(g.location_name == loc) & (g.year == yr) & g.age_name.eq('20+ years') & g.sex_name.eq('Both')]
        res[f'gbd_{loc}_{yr}_20plus'] = {m: [float(d[d.measure_name == m].val.iloc[0]), float(d[d.measure_name == m].lower.iloc[0]), float(d[d.measure_name == m].upper.iloc[0])] for m in ['Prevalence', 'Deaths']}
json.dump(res, open(out, 'w'), indent=1); print(json.dumps(res, indent=1))
