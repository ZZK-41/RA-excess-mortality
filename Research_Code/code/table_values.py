"""Values for the Canadian table, the European comparisons (table S7), and remaining text numbers. Args: ccdss_raw statcan_zip atlas_dir who_csv t5_csv out_json"""
import sys, io, json, zipfile
import numpy as np, pandas as pd
ccdss, statcan, atlas, whocsv, t5p, out = sys.argv[1:7]; Z = 1.959964; B = 20000; SEED = 20261123
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
for a, b in [('d1', 'Counts of mortality with the condition'), ('n1', 'Population of mortality with the condition'), ('d0', 'Counts of mortality without the condition'), ('n0', 'Population of mortality without the condition')]: t[a] = nm(t[b])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[c5.age >= 40].copy()
c5['sex'] = c5.Sex.map({'Females': 'Women', 'Males': 'Men'}); c5['N'] = c5.n1 + c5.n0; c5['mpop'] = (c5.d1 + c5.d0)/c5.N; c5['E'] = c5.n1*c5.mpop
w = c5[c5.year.between(2007, 2009)].groupby(['age', 'sex']).n1.sum(); w = (w/w.sum()).rename('w').reset_index(); c5 = c5.merge(w, on=['age', 'sex'])
sc = pd.read_csv(io.BytesIO(zipfile.ZipFile(statcan).read('13100150.csv')), low_memory=False)
sc = sc[sc['Cause of death (ICD-10)'].isin(['Seropositive rheumatoid arthritis [M05]', 'Other rheumatoid arthritis [M06]'])]
amap = {f'{a} to {a+4} years': a for a in range(40, 90, 5)} | {'90 years and over': 90}
sc['age'] = sc['Age group'].map(amap); sc['sex'] = sc.Sex.map({'Females': 'Women', 'Males': 'Men', 'Both sexes': 'Both'}); sc = sc[sc.age.notna()]
cert = sc[sc.sex.isin(['Women', 'Men'])].groupby(['REF_DATE', 'age', 'sex']).VALUE.sum().rename('C').reset_index().rename(columns={'REF_DATE': 'year'})
c5 = c5.merge(cert, on=['year', 'age', 'sex'], how='left').fillna({'C': 0})
assert abs(c5[c5.year == 2023].C.sum() - 331) < 0.5
def row(d):
    D, E, C, n = d.d1.sum(), d.E.sum(), d.C.sum(), d.n1.sum(); rng = np.random.default_rng(SEED)
    r = (rng.poisson(D, B) - E)/np.maximum(rng.poisson(C, B), 1)
    return dict(people=float(n), deaths=float(D), expected=float(E), excess=float(D - E), excess_lo=float(D - E - Z*np.sqrt(D)), excess_hi=float(D - E + Z*np.sqrt(D)),
                smr=float(D/E), smr_lo=float((D - Z*np.sqrt(D))/E), smr_hi=float((D + Z*np.sqrt(D))/E), attributed=float(C), ratio=float((D - E)/C),
                ratio_lo=float(np.quantile(r, .025)), ratio_hi=float(np.quantile(r, .975)), attributed_pct_of_deaths=float(100*C/D), attributed_pct_of_excess=float(100*C/(D - E)))
K = {'table1': {}}
for per, yrs in [('2023', [2023]), ('2007-23', list(range(2007, 2024)))]:
    d = c5[c5.year.isin(yrs)]
    K['table1'][per] = {'All ≥40': row(d), '40–64': row(d[d.age < 65]), '65–79': row(d[(d.age >= 65) & (d.age < 80)]), '≥80': row(d[d.age >= 80]),
                        'Women': row(d[d.sex == 'Women']), 'Men': row(d[d.sex == 'Men'])}
# standardised rates: % change 2007-09 -> 2021-23 with simulation CIs
rng = np.random.default_rng(SEED); Bs = 5000
def std(g, D1, C):
    return (1000*(g.w.values*C/g.n1.values)).sum(), (1000*(g.w.values*(D1/g.n1.values - (D1 + g.d0.values)/g.N.values))).sum()
res = {'cert': [], 'exc': [], 'share': []}
p1 = c5[c5.year.between(2007, 2009)]; p2 = c5[c5.year.between(2021, 2023)]
pt = {}
for lab, p in [('p1', p1), ('p2', p2)]:
    vals = [std(g, g.d1.values, g.C.values) for _, g in p.groupby('year')]; pt[lab] = np.mean(vals, axis=0)
for b in range(Bs):
    v = {}
    for lab, p in [('p1', p1), ('p2', p2)]:
        vals = [std(g, rng.poisson(g.d1.values), rng.poisson(g.C.values)) for _, g in p.groupby('year')]; v[lab] = np.mean(vals, axis=0)
    res['cert'].append(100*(v['p2'][0]/v['p1'][0] - 1)); res['exc'].append(100*(v['p2'][1]/v['p1'][1] - 1)); res['share'].append(100*((v['p2'][0]/v['p2'][1])/(v['p1'][0]/v['p1'][1]) - 1))
K['std_change'] = dict(certified=[float(100*(pt['p2'][0]/pt['p1'][0] - 1))] + [float(x) for x in np.quantile(res['cert'], [.025, .975])],
                       excess=[float(100*(pt['p2'][1]/pt['p1'][1] - 1))] + [float(x) for x in np.quantile(res['exc'], [.025, .975])],
                       share=[float(100*((pt['p2'][0]/pt['p2'][1])/(pt['p1'][0]/pt['p1'][1]) - 1))] + [float(x) for x in np.quantile(res['share'], [.025, .975])],
                       levels=dict(cert_0709=float(pt['p1'][0]), cert_2123=float(pt['p2'][0]), exc_0709=float(pt['p1'][1]), exc_2123=float(pt['p2'][1])))
# Denmark Atlas: deaths, people, MRR with CI; registered deaths 2000-18
me = pd.read_csv(f'{atlas}/atlas_02_main_estimates.txt', sep=' ', quotechar='"'); mr = pd.read_csv(f'{atlas}/atlas_04_MRR.txt', sep=' ', quotechar='"')
who = pd.read_csv(whocsv); dk = who[(who.country == 'Denmark') & who.Year.between(2000, 2018)]
K['denmark'] = {'registered_2000_18': float(dk.M05_M06.sum()), 'registered_M05': float(dk.M05.sum()), 'registered_M06': float(dk.M06.sum()), 'years': int(len(dk))}
for code, lab in [(10919, 'M05'), (10920, 'M06')]:
    m = me[(me.id == code) & (me.sex == 'All')].iloc[0]; r = mr[(mr.id == code) & (mr.sex == 'All') & (mr.model == 'All')].iloc[0]
    m = m.copy(); m['n_deaths'] = float(str(m.n_deaths).replace(',', '')); m['prev_cases'] = float(str(m.prev_cases).replace(',', ''))
    r = r.copy(); r['MRR'], r['MRR_left'], r['MRR_right'] = float(r.MRR), float(r.MRR_left), float(r.MRR_right)
    K['denmark'][lab] = dict(people=float(m.prev_cases), deaths=float(m.n_deaths), mrr=float(r.MRR), mrr_lo=float(r.MRR_left), mrr_hi=float(r.MRR_right),
                             excess=float(m.n_deaths*(1 - 1/r.MRR)), excess_lo=float(m.n_deaths*(1 - 1/r.MRR_left)), excess_hi=float(m.n_deaths*(1 - 1/r.MRR_right)))
    print(lab, {k: v for k, v in m.items() if k in ('n_deaths', 'prev_cases', 'age_death_median', 'age_diag_median')})
K['france_registered_2019'] = float(who[(who.country == 'France') & (who.Year == 2019)].M05_M06.iloc[0])
t5 = pd.read_csv(t5p); fr = t5[t5.setting.str.contains('France')].iloc[0]; K['france'] = {k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in fr.items()}
json.dump(K, open(out, 'w'), indent=1, ensure_ascii=False)
for per in K['table1']:
    for k, v in K['table1'][per].items(): print(per, k, {kk: round(vv, 3) for kk, vv in v.items()})
print(json.dumps({k: K[k] for k in ['std_change', 'denmark', 'france_registered_2019']}, indent=1))
