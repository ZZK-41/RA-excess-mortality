"""Build the GBD 2023 country cells with the four mortality profiles (ages 45-89 years, 2019 and 2023).

Inputs: the three GBD 2023 results-tool files for 2019 and 2023 (RA and all-cause measures for 204 countries and the world),
the CCDSS export (Canadian rate ratios and rate differences for people with and without RA), and the Danish atlas file of
age-specific mortality (atlas_06_incidences_by_age.txt; M05 = id 10919, M06 = id 10920, whole population = id 0).
For every location, year, sex, and 5-year age group the mortality of people with RA (m1) is derived from GBD all-cause mortality (m)
and RA prevalence (p) under each profile, and modelled excess deaths are GBD cases x (m1 - m).
The GBD values for all ages are first written to an intermediate file (gbd_analysis_cells.csv) and read back, as in the original
analysis, so that results are reproduced to the last digit.
Usage: python3 -I build_gbd_cells.py <gbd_2019_2023_dir> <ccdss_csv> <atlas_06_txt> <out_csv_gz> <intermediate_csv>
"""
import sys, io, zipfile, itertools
from pathlib import Path
import numpy as np, pandas as pd

gdir, ccdss, atlas06, out, inter = sys.argv[1:6]
FILES = sorted(['IHME-GBD_2023_DATA-aaf6002d-1.zip', 'IHME-GBD_2023_DATA-c82898aa-1.zip', 'IHME-GBD_2023_DATA-8e54c028-1.zip'])
frames = []
for f in FILES:
    z = zipfile.ZipFile(Path(gdir) / f)
    for n in z.namelist():
        if n.endswith('.csv'):
            frames.append(pd.read_csv(io.BytesIO(z.read(n))))
g = pd.concat(frames, ignore_index=True)
keys = ['cause_id', 'measure_id', 'metric_id', 'location_id', 'age_id', 'sex_id', 'year']
assert not g.duplicated(keys).any()
ages = list(range(9, 21)) + [30, 31, 32, 235, 37]
locs = sorted(g.location_id.unique()); assert len(locs) == 205 and 1 in locs
need = set(itertools.product([627], [1, 3, 4, 5], [1, 3], locs, ages, [1, 2, 3], [2019, 2023])) | set(itertools.product([294], [1], [1, 3], locs, ages, [1, 2, 3], [2019, 2023]))
have = set(g[keys].itertuples(index=False, name=None)); assert not need - have, 'missing GBD cells'
idx = ['location_id', 'location_name', 'year', 'sex_id', 'sex_name', 'age_id', 'age_name']
wide = g.pivot(index=idx, columns=['cause_id', 'measure_id', 'metric_id'], values='val')
a = pd.DataFrame(index=wide.index)
a['cases'] = wide[(627, 5, 1)]; a['p'] = wide[(627, 5, 3)] / 1e5; a['m'] = wide[(294, 1, 3)] / 1e5
a['pop_ra'] = a.cases / a.p; a['pop_deaths'] = wide[(294, 1, 1)] / a.m
a['ra_deaths'] = wide[(627, 1, 1)]; a['ylds'] = wide[(627, 3, 1)]; a['ylls'] = wide[(627, 4, 1)]
ra = g[(g.cause_id == 627) & (g.measure_id == 5) & (g.metric_id == 1)].set_index(idx)
a['cases_lower'] = ra.lower; a['cases_upper'] = ra.upper
bg = g[(g.cause_id == 294) & (g.measure_id == 1) & (g.metric_id == 3)].set_index(idx)
a['m_lower'] = bg.lower / 1e5; a['m_upper'] = bg.upper / 1e5
a = a.reset_index(); a['age_start'] = a.age_name.str.extract(r'^(\d+)')[0].astype(float)
a.to_csv(inter, index=False); a = pd.read_csv(inter)
a = a[a.sex_id.isin([1, 2]) & a.age_start.between(20, 95) & ~a.age_id.isin([37, 22])]
x = a[a.age_start.between(45, 85)].copy().sort_values(['location_id', 'year', 'sex_id', 'age_start'])
# Canadian profiles: CCDSS national 5-year cells, fiscal years 2019-20 and 2023-24
raw = open(ccdss, encoding='latin-1').read().splitlines()
c = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); c.columns = [k.replace('\xa0', ' ') for k in c.columns]
c = c[c['Fiscal year'].notna() & c.Sex.notna()].copy()
for k in ['Age group', 'Fiscal year']: c[k] = c[k].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
c = c[(c.Geography == 'Canada') & (c['Age group type'] == '5 Year') & c.Sex.isin(['Females', 'Males'])].copy()
num = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
c['d1'] = num(c['Counts of mortality with the condition']); c['n1'] = num(c['Population of mortality with the condition'])
c['d0'] = num(c['Counts of mortality without the condition']); c['n0'] = num(c['Population of mortality without the condition'])
c['year'] = c['Fiscal year'].str[:4].astype(int); c['sex_id'] = c.Sex.map({'Males': 1, 'Females': 2})
c['age_start'] = pd.to_numeric(c['Age group'].str.extract(r'^(\d+)')[0])
c = c[c.year.isin([2019, 2023]) & c.age_start.between(45, 85)].copy()
c['R'] = (c.d1 / c.n1) / (c.d0 / c.n0); c['delta0'] = c.d1 / c.n1 - c.d0 / c.n0
c['ratio_gp'] = (c.d1 / c.n1) / ((c.d1 + c.d0) / (c.n1 + c.n0))
assert len(c) == 36 and c.R.notna().all()
x = x.merge(c[['year', 'sex_id', 'age_start', 'R', 'delta0', 'ratio_gp', 'd1', 'd0', 'n1', 'n0']], on=['year', 'sex_id', 'age_start'], validate='many_to_one')
# Danish profiles: mortality of people with a hospital diagnosis of M05 or M06 relative to the whole population, by sex and age
t = pd.read_csv(atlas06, sep=' ', quotechar='"', low_memory=False)
t = t[t.id.isin([0, 10919, 10920])].copy()
t['age_start'] = pd.to_numeric(t.age_cat.str.extract(r'\[(\d+)')[0]); t['sex_id'] = t.sex.map({'Males': 1, 'Females': 2})
for ident, col in [(10919, 'Danish M05'), (10920, 'Danish M06')]:
    rr = t[t.id == ident][['sex_id', 'age_start', 'death_rate']].rename(columns={'death_rate': 'r'}).merge(
        t[t.id == 0][['sex_id', 'age_start', 'death_rate']].rename(columns={'death_rate': 'b'}), on=['sex_id', 'age_start'])
    rr[col] = rr.r / rr.b
    x = x.merge(rr[['sex_id', 'age_start', col]], on=['sex_id', 'age_start'], validate='many_to_one')
assert x[['R', 'Danish M05', 'Danish M06', 'm', 'p']].notna().all().all()
parts = []
for scenario in ['Canada relative', 'Canada additive', 'Danish M05', 'Danish M06', 'Canada 2019 fixed']:
    z = x.copy(); z['scenario'] = scenario
    if scenario == 'Canada relative': z['m1'] = z.m * z.R / (1 + z.p * (z.R - 1))
    elif scenario == 'Canada additive': z['m1'] = z.m + (1 - z.p) * z.delta0
    elif scenario.startswith('Danish'): z['m1'] = z.m * z[scenario]
    else:
        v = c[c.year == 2019][['sex_id', 'age_start', 'R']].rename(columns={'R': 'R_fixed'})
        z = z.merge(v, on=['sex_id', 'age_start'], validate='many_to_one'); z['m1'] = z.m * z.R_fixed / (1 + z.p * (z.R_fixed - 1))
    z['difference'] = z.m1 - z.m; z['equivalent'] = z.cases * z.difference; z['m0_implied'] = (z.m - z.p * z.m1) / (1 - z.p)
    assert (z.m1 > 0).all() and (z.m0_implied > 0).all()
    parts.append(z)
cells = pd.concat(parts, ignore_index=True)
cells.to_csv(out, index=False)
print('GBD cells:', len(cells), 'rows,', cells.location_id.nunique(), 'locations')
