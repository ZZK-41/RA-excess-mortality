"""Remaining numbers: additional deaths in 2020-23 above the FY2015-19 relative excess (ages >=40),
non-RA comparator sensitivity, age shares, per-person excess trend. Args: ccdss_raw out_json"""
import sys, io, json
import numpy as np, pandas as pd
raw = open(sys.argv[1], encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
t['d1'] = nm(t['Counts of mortality with the condition']); t['n1'] = nm(t['Population of mortality with the condition'])
t['d0'] = nm(t['Counts of mortality without the condition']); t['n0'] = nm(t['Population of mortality without the condition'])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[c5.age >= 40]
c5['N'] = c5.n1 + c5.n0; c5['mp'] = (c5.d1 + c5.d0)/c5.N; c5['E'] = c5.n1*c5.mp; c5['E0'] = c5.n1*c5.d0/c5.n0
base = c5[c5.year.between(2015, 2019)].groupby(['age', 'Sex']).apply(lambda g: (g.d1.sum()/g.n1.sum())/((g.d1+g.d0).sum()/g.N.sum()), include_groups=False).rename('rr').reset_index()
out = {}
Z = 1.959964
for yr in [2020, 2021, 2022, 2023]:
    g = c5[c5.year == yr].merge(base, on=['age', 'Sex']); e = (g.n1*g.mp*g.rr).sum(); o = g.d1.sum()
    out[f'additional_{yr}'] = [o - e, o - e - Z*np.sqrt(o), o - e + Z*np.sqrt(o)]
s = c5[c5.year == 2023]
out['excess_vs_nonra_2023'] = s.d1.sum() - s.E0.sum(); out['excess_vs_pop_2023'] = s.d1.sum() - s.E.sum()
out['share_deaths_ge80_2023'] = s[s.age >= 80].d1.sum()/s.d1.sum(); out['share_excess_ge80_2023'] = (s[s.age >= 80].d1.sum() - s[s.age >= 80].E.sum())/(s.d1.sum() - s.E.sum())
out['prevalence_ge40_2023_pct'] = 100*s.n1.sum()/s.N.sum(); out['people_ge80_share_2023'] = s[s.age >= 80].n1.sum()/s.n1.sum()
a = c5.groupby('year').apply(lambda g: pd.Series(dict(n=g.n1.sum(), ex=(g.d1 - g.E).sum(), per1000=1000*(g.d1 - g.E).sum()/g.n1.sum())), include_groups=False)
out['per1000_2007'] = a.loc[2007, 'per1000']; out['per1000_2019'] = a.loc[2019, 'per1000']; out['per1000_2023'] = a.loc[2023, 'per1000']
out['people_2007'] = a.loc[2007, 'n']; out['people_2023'] = a.loc[2023, 'n']; out['excess_2007'] = a.loc[2007, 'ex']
out['people_ge80_2007'] = c5[(c5.year == 2007) & (c5.age >= 80)].n1.sum(); out['people_ge80_2023'] = s[s.age >= 80].n1.sum()
json.dump(out, open(sys.argv[2], 'w'), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
