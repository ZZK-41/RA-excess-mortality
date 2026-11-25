"""Robustness of the world observed-vs-expected comparison (see gbd_expected.py) to the standard population, age range, and end year.
Units: 7 GBD super-regions. Args: gbd_annual_dir out_json"""
import sys, io, json, zipfile
import pandas as pd
gbd_annual, out = sys.argv[1:3]
def rz(p):
    z = zipfile.ZipFile(p); return pd.concat([pd.read_csv(io.BytesIO(z.read(n)), low_memory=False) for n in z.namelist() if n.endswith('.csv')])
AGE40 = [13, 14, 15, 16, 17, 18, 19, 20, 30, 31, 32, 235]; SR = [4, 31, 64, 103, 137, 158, 166]; key = ['location_id', 'sex_id', 'age_id', 'year']; cell = ['sex_id', 'age_id']
ra = pd.concat([rz(f'{gbd_annual}/IHME-GBD_2023_DATA-{x}-1.zip') for x in ['04c719ed', 'd1117f20', 'f5ad253b']]).drop_duplicates(['measure_id', 'metric_id', 'location_id', 'sex_id', 'age_id', 'year'])
pv = pd.concat([rz(f'{gbd_annual}/IHME-GBD_2023_DATA-d9e5f203-1.zip'), rz(f'{gbd_annual}/IHME-GBD_2023_DATA-a4d3c489-1.zip')]); pv = pv[pv.measure_id == 5]
ac = rz(f'{gbd_annual}/IHME-GBD_2023_DATA-0e1e6620-1.zip')
D = ra[(ra.measure_id == 1) & (ra.metric_id == 1) & ra.sex_id.isin([1, 2]) & ra.age_id.isin(AGE40)].set_index(key).val.rename('d')
C = pv[(pv.metric_id == 1) & pv.sex_id.isin([1, 2]) & pv.age_id.isin(AGE40)].set_index(key).val.rename('c')
M = (ac[ac.sex_id.isin([1, 2]) & ac.age_id.isin(AGE40) & (ac.metric_id == 3)].set_index(key).val/1e5).rename('m')
an = pd.concat([D, C, M], axis=1).dropna()
def run(wyear, ages, t0=1990, t1=2023):
    x = an[an.index.get_level_values('age_id').isin(ages) & an.index.get_level_values('location_id').isin(SR)]
    w = an.xs((1, wyear), level=['location_id', 'year']).c; w = w[w.index.get_level_values('age_id').isin(ages)]; w = w/w.sum()
    b = x.xs(t0, level='year'); a0 = b.d/b.c
    def at(t):
        g = x.xs(t, level='year'); pi = g.c/g.c.groupby(level=cell).transform('sum')
        o = (pi*g.d/g.c).groupby(level=cell).sum(); e = (pi*a0.reindex(g.index)*g.m/b.m.reindex(g.index)).groupby(level=cell).sum()
        return (w*o).sum(), (w*e).sum()
    o0, _ = at(t0); o1, e1 = at(t1)
    return dict(observed=float(100*(o1/o0 - 1)), expected=float(100*(e1/o0 - 1)), observed_vs_expected=float(100*(o1/e1 - 1)))
R = {'main_2019_weights_ge40_1990_2023': run(2019, AGE40), 'weights_1990': run(1990, AGE40), 'weights_2023': run(2023, AGE40),
     'ages_40_79': run(2019, AGE40[:8]), 'ages_40_89': run(2019, AGE40[:10]), 'end_2019': run(2019, AGE40, t1=2019), 'start_2000_end_2019': run(2019, AGE40, t0=2000, t1=2019)}
json.dump(R, open(out, 'w'), indent=1)
for k, v in R.items(): print(f"{k:36s} observed {v['observed']:6.1f}  expected {v['expected']:6.1f}  observed vs expected {v['observed_vs_expected']:6.1f}")
