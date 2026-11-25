"""Registered underlying-cause RA deaths (ICD-10 M05, M06) from WHO Mortality Database raw files.
Args: out_csv zip1 [zip2 ...]. Reads data only; no code from the archives is executed."""
import sys, io, zipfile
import pandas as pd
out = sys.argv[1]; zips = sys.argv[2:]
KEEP = {2090: 'Canada', 4050: 'Denmark', 4080: 'France', 4180: 'Italy', 4290: 'Sweden', 4220: 'Norway', 4070: 'Finland', 4210: 'Netherlands', 4308: 'United Kingdom', 2450: 'United States of America', 5020: 'Australia', 3160: 'Japan'}
rows = []
for p in zips:
    z = zipfile.ZipFile(p)
    for n in z.namelist():
        for ch in pd.read_csv(io.BytesIO(z.read(n)), dtype={'Cause': str, 'List': str}, usecols=['Country', 'Admin1', 'SubDiv', 'Year', 'List', 'Cause', 'Sex', 'Frmat', 'Deaths1'], chunksize=500000, low_memory=False):
            ch = ch[ch.Country.isin(KEEP) & ch.Admin1.isna() & ch.SubDiv.isna()]
            c = ch.Cause.astype(str).str.upper()
            ra = ch[c.str.startswith('M05') | c.str.startswith('M06')].copy()
            allc = ch[c.isin(['AAA', '1000'])].copy()
            if len(ra): ra['code'] = ra.Cause.str[:3]; rows.append(ra.assign(kind='ra'))
            if len(allc): allc['code'] = 'ALL'; rows.append(allc.assign(kind='all'))
d = pd.concat(rows); d['country'] = d.Country.map(KEEP)
# keep one tabulation list per country-year for RA: prefer detailed lists (104, 103, 10M); 'ALL' from any list
ra = d[d.kind == 'ra']; best = ra.groupby(['country', 'Year']).List.agg(lambda s: sorted(set(s))[0]); 
ra = ra.merge(best.rename('L').reset_index(), on=['country', 'Year']); ra = ra[ra.List == ra.L]
s = ra[ra.Sex.isin([1, 2])].groupby(['country', 'Year', 'List', 'code']).Deaths1.sum().unstack('code').fillna(0).reset_index()
s['M05_M06'] = s.get('M05', 0) + s.get('M06', 0)
al = d[(d.kind == 'all') & d.Sex.isin([1, 2])].groupby(['country', 'Year', 'List']).Deaths1.sum().groupby(['country', 'Year']).max().rename('all_cause').reset_index()
s = s.merge(al, on=['country', 'Year'], how='left'); s.to_csv(out, index=False)
pd.set_option('display.width', 200); pd.set_option('display.max_rows', 400)
print(s[s.country.isin(['Canada', 'Denmark', 'France', 'Italy'])].to_string(index=False))
print(s.groupby('country').Year.agg(['min', 'max', 'count']))
