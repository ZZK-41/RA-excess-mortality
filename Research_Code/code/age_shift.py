"""Attributed share of excess deaths by age band, 2007-09 vs 2021-23, and age-standardised counterfactual.
Args: ccdss_raw statcan_zip [out_csv]"""
import sys, io, zipfile
import numpy as np, pandas as pd
raw = open(sys.argv[1], encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
for a, b in [('d1', 'Counts of mortality with the condition'), ('n1', 'Population of mortality with the condition'), ('d0', 'Counts of mortality without the condition'), ('n0', 'Population of mortality without the condition')]: t[a] = nm(t[b])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[c5.age >= 40]
c5['ex'] = c5.d1 - c5.n1*(c5.d1 + c5.d0)/(c5.n1 + c5.n0)
sc = pd.read_csv(io.BytesIO(zipfile.ZipFile(sys.argv[2]).read('13100150.csv')), low_memory=False)
sc = sc[sc['Cause of death (ICD-10)'].isin(['Seropositive rheumatoid arthritis [M05]', 'Other rheumatoid arthritis [M06]']) & sc.Sex.eq('Both sexes')]
amap = {f'{a} to {a+4} years': a for a in range(40, 90, 5)} | {'90 years and over': 90}
sc['age'] = sc['Age group'].map(amap); sc = sc[sc.age.notna()]
bands = {'40–64': (40, 65), '65–79': (65, 80), '80–89': (80, 90), '≥90': (90, 200)}
rows = []
for lab, yrs in [('2007–09', [2007, 2008, 2009]), ('2021–23', [2021, 2022, 2023])]:
    for b, (lo, hi) in bands.items():
        ex = c5[c5.year.isin(yrs) & (c5.age >= lo) & (c5.age < hi)].ex.sum(); co = sc[sc.REF_DATE.isin(yrs) & (sc.age >= lo) & (sc.age < hi)].VALUE.sum()
        rows.append(dict(period=lab, ages=b, excess=ex, coded=co, coded_pct=100*co/ex))
r = pd.DataFrame(rows); print(r.round(2).to_string(index=False))
if len(sys.argv) > 3: r.to_csv(sys.argv[3], index=False)
p = r.pivot(index='ages', columns='period', values='excess'); s = r.pivot(index='ages', columns='period', values='coded_pct')/100
w07 = p['2007–09']/p['2007–09'].sum(); w21 = p['2021–23']/p['2021–23'].sum()
print('overall coded share 2007-09 %.2f%%, 2021-23 %.2f%%' % (100*r[r.period == '2007–09'].coded.sum()/r[r.period == '2007–09'].excess.sum(), 100*r[r.period == '2021–23'].coded.sum()/r[r.period == '2021–23'].excess.sum()))
print('2021-23 share with 2007-09 age distribution of excess: %.2f%%' % (100*(w07*s['2021–23']).sum()))
print('2007-09 share with 2021-23 age distribution of excess: %.2f%%' % (100*(w21*s['2007–09']).sum()))
print('excess age distribution 2007-09:', w07.round(3).to_dict(), ' 2021-23:', w21.round(3).to_dict())
