"""Checks of the GBD decomposition and of annual tracking of certified deaths. Args: gbd_annual_dir gbd_2019_2023_dir statcan_zip trends_dir"""
import sys, io, zipfile, json
import numpy as np, pandas as pd
gbd_annual, gbd_2019_2023, statcan, td = sys.argv[1:5]
def rz(p):
    z = zipfile.ZipFile(p); return pd.concat([pd.read_csv(io.BytesIO(z.read(n)), low_memory=False) for n in z.namelist() if n.endswith('.csv')])
P2F, P2G, P2H = [rz(f'{gbd_annual}/IHME-GBD_2023_DATA-{x}-1.zip') for x in ['80a09abe', '4af37d63', '532d5fea']]
R19, R23, A = [rz(f'{gbd_2019_2023}/IHME-GBD_2023_DATA-{x}-1.zip') for x in ['c82898aa', 'aaf6002d', '8e54c028']]
ann = pd.concat([rz(f'{gbd_annual}/IHME-GBD_2023_DATA-{x}-1.zip') for x in ['04c719ed', 'd1117f20', 'f5ad253b']]).drop_duplicates(['measure_id', 'metric_id', 'location_id', 'sex_id', 'age_id', 'year'])
def cells(yr, ages):
    ra = {2000: P2F, 2010: P2G, 2019: R19, 2023: R23}[yr]; ac = P2H if yr in (2000, 2010) else A
    ra = ra[(ra.metric_id == 1) & ra.sex_id.isin([1, 2]) & ra.age_id.isin(ages) & (ra.year == yr)]
    p = ra.pivot_table(index=['location_id', 'sex_id', 'age_id'], columns='measure_id', values='val').rename(columns={1: 'd', 5: 'c'})
    a = ac[ac.sex_id.isin([1, 2]) & ac.age_id.isin(ages) & (ac.year == yr)].pivot_table(index=['location_id', 'sex_id', 'age_id'], columns='metric_id', values='val')
    p['m'] = a[3]/1e5; return p.reset_index()
def decomp(loc, ages, wyear):
    cs = {y: cells(y, ages) for y in [2000, 2010, 2019, 2023]}
    w = cs[wyear][cs[wyear].location_id == 1].set_index(['sex_id', 'age_id']).c; w = w/w.sum()
    out = {}
    for y, c in cs.items():
        if isinstance(loc, int): q = c[c.location_id == loc].set_index(['sex_id', 'age_id'])
        else: q = c[c.location_id.isin(loc)].groupby(['sex_id', 'age_id']).apply(lambda g: pd.Series(dict(d=g.d.sum(), c=g.c.sum(), m=(g.m*g.c).sum()/g.c.sum())), include_groups=False)
        out[y] = ((w*q.d/q.c).sum(), (w*q.m).sum())
    att = 100*(out[2023][0]/out[2000][0] - 1); bg = 100*(out[2023][1]/out[2000][1] - 1)
    att10 = 100*(out[2023][0]/out[2010][0] - 1); bg10 = 100*(out[2023][1]/out[2010][1] - 1)
    return f'RA deaths per case {att:+.1f}%, all-cause {bg:+.1f}% (2000-23); 2010-23: {att10:+.1f}% vs {bg10:+.1f}%'
A40 = [13, 14, 15, 16, 17, 18, 19, 20, 30, 31, 32, 235]; A20 = [9, 10, 11, 12] + A40
ages_avail = sorted(set(P2F.age_id))
print('ages in 2000 file', ages_avail)
for lab, loc in [('World', 1), ('High-income', 64), ('Canada', 101)]:
    print(lab, '| ages>=40, 2019 weights:', decomp(loc, A40, 2019)); print(lab, '| ages>=40, 2000 weights:', decomp(loc, A40, 2000))
    print(lab, '| ages 40-79:', decomp(loc, [13, 14, 15, 16, 17, 18, 19, 20], 2019)); print(lab, '| ages >=80:', decomp(loc, [30, 31, 32, 235], 2019))
# Canada annual GBD vs StatCan certified, ages >=40
g = ann[(ann.measure_id == 1) & (ann.metric_id == 1) & (ann.location_id == 101) & (ann.sex_id == 3) & ann.age_id.isin(A40)].groupby('year').val.sum()
sc = pd.read_csv(io.BytesIO(zipfile.ZipFile(statcan).read('13100150.csv')), low_memory=False)
sc = sc[sc['Cause of death (ICD-10)'].isin(['Seropositive rheumatoid arthritis [M05]', 'Other rheumatoid arthritis [M06]']) & sc.Sex.eq('Both sexes')]
sc['age'] = sc['Age group'].map({f'{x} to {x+4} years': x for x in range(40, 90, 5)} | {'90 years and over': 90}); c = sc[sc.age.notna()].groupby('REF_DATE').VALUE.sum()
d = pd.DataFrame({'gbd': g, 'certified': c}).dropna(); d = d.loc[2000:2023]
print('Canada GBD vs certified >=40, 2000-23: ratio median %.3f (%.3f-%.3f); r levels %.3f; r annual changes %.3f' % ((d.gbd/d.certified).median(), (d.gbd/d.certified).min(), (d.gbd/d.certified).max(), np.corrcoef(d.gbd, d.certified)[0, 1], np.corrcoef(np.diff(d.gbd), np.diff(d.certified))[0, 1]))
print(d.loc[2019:2023].round(1).to_string())
# YLL share of GBD RA DALYs, global, ages >=20
y = ann[(ann.metric_id == 1) & (ann.location_id == 1) & (ann.sex_id == 3) & (ann.age_id == 37)].pivot_table(index='year', columns='measure_id', values='val')
for yr in [1990, 2000, 2010, 2023]: print(yr, 'YLL share of RA DALYs (global, 20+): %.1f%%' % (100*y.loc[yr, 4]/y.loc[yr, 2]))
