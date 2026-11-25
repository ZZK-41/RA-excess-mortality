"""Build the provincial analysis panel from the CCDSS data-tool export.

Rows: province x fiscal year (2010-11 onwards) x age group (50-64, 65-79, >=80 years) x sex, kept only when all six
age-sex strata of the province-year are published with positive counts. The weight column is the number of people
with RA in the same stratum, summed over the province-years 2007-08 to 2009-10 that meet the same rule.
Usage: python3 -I build_provincial_panel.py <ccdss_csv> <out_csv>
"""
import sys, io
import pandas as pd

src, out = sys.argv[1:3]
raw = open(src, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False)          # first line is the table title
t.columns = [c.replace('\xa0', ' ') for c in t.columns]
t = t[t['Fiscal year'].notna() & t.Sex.notna()].copy()                       # drop the notes below the table
for c in ['Age group', 'Fiscal year']:
    t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
num = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
COLS = {'d_ra': 'Counts of mortality with the condition', 'n_ra': 'Population of mortality with the condition',
        'd_nonra': 'Counts of mortality without the condition', 'n_nonra': 'Population of mortality without the condition'}
for k, c in COLS.items():
    t[k] = num(t[c])
t['year'] = t['Fiscal year'].str[:4].astype(int)
p = t[(t.Geography != 'Canada') & t.Sex.isin(['Females', 'Males']) & t['Age group'].isin(['50–64', '65–79', '80+'])].copy()
p['published'] = p[list(COLS)].notna().all(axis=1) & (p[list(COLS)] > 0).all(axis=1)
p['complete'] = p.groupby(['Geography', 'year']).published.transform('sum').eq(6)
p = p[p.complete].rename(columns={'Geography': 'geography', 'Age group': 'age_group', 'Sex': 'sex'})
w = p[p.year < 2010].groupby(['age_group', 'sex']).n_ra.sum().rename('weight')
panel = p[p.year >= 2010].merge(w.reset_index(), on=['age_group', 'sex'], how='left')
panel['older'] = panel.age_group.eq('80+')
panel = panel[['geography', 'year', 'age_group', 'sex', 'd_ra', 'd_nonra', 'n_ra', 'n_nonra', 'weight', 'older']].sort_values(['geography', 'year', 'age_group', 'sex'])
for k in ['d_ra', 'd_nonra', 'n_ra', 'n_nonra', 'weight']:
    panel[k] = panel[k].astype(int)
assert len(panel) == 540 and panel.weight.notna().all()
panel.to_csv(out, index=False)
print('provincial panel:', len(panel), 'rows;', panel.groupby('geography').year.nunique().to_dict())
