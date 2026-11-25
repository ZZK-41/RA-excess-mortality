"""Canada 2023, ages >=40: GBD RA deaths/YLL/YLD versus observed excess deaths and their YLL (CCDSS case numbers and GBD case numbers).
YLL per death from GBD 2023 Canada by sex and 5-year age (ages 90+ pooled). 95% CIs: Poisson on observed RA deaths.
Args: ccdss_raw gbd_dir out_json"""
import sys, io, json, zipfile
import numpy as np, pandas as pd
ccdss, gbddir, out = sys.argv[1:4]; Z = 1.959964; B = 20000
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False); t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
for a, b in [('d1', 'Counts of mortality with the condition'), ('n1', 'Population of mortality with the condition'), ('d0', 'Counts of mortality without the condition'), ('n0', 'Population of mortality without the condition')]: t[a] = nm(t[b])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int); c5 = c5[(c5.age >= 40) & (c5.year == 2023)].copy()
c5['sex'] = c5.Sex.map({'Females': 'Female', 'Males': 'Male'}); c5['E'] = c5.n1*(c5.d1 + c5.d0)/(c5.n1 + c5.n0)
z = zipfile.ZipFile(f'{gbddir}/IHME-GBD_2023_DATA-aaf6002d-1.zip'); m = [m for m in z.infolist() if m.filename.endswith('.csv')][0]
g = pd.read_csv(io.BytesIO(z.read(m)), low_memory=False)
g = g[(g.location_name == 'Canada') & (g.metric_name == 'Number') & g.sex_name.isin(['Female', 'Male'])]
FIVE = {f'{a}-{a+4} years': a for a in range(20, 95, 5)} | {'95+ years': 95}
g = g[g.age_name.isin(FIVE)].copy(); g['age'] = g.age_name.map(FIVE).clip(upper=90)
p = g.pivot_table(index=['sex_name', 'age'], columns='measure_name', values='val', aggfunc='sum').reset_index().rename(columns={'sex_name': 'sex'})
p = p.rename(columns={'YLLs (Years of Life Lost)': 'yll', 'YLDs (Years Lived with Disability)': 'yld', 'Deaths': 'deaths', 'Prevalence': 'cases'}); p['ypd'] = p.yll/p.deaths
s = c5.merge(p, on=['sex', 'age']); assert len(s) == 22
f = s.cases/s.n1
rng = np.random.default_rng(20261123); D = rng.poisson(s.d1.values, (B, len(s)))
ex_sim = D - s.E.values; yll_full = (ex_sim*s.ypd.values).sum(1); yll_cons = (ex_sim*(f*s.ypd).values).sum(1)
gbd40 = dict(deaths=float(s.deaths.sum()), yll=float(s.yll.sum()), yld=float(s.yld.sum()), cases=float(s.cases.sum()))
gbd40['daly'] = gbd40['yll'] + gbd40['yld']
res = dict(gbd_ge40=gbd40)
for lab, sim, w in [('ccdss', yll_full, np.ones(len(s))), ('gbd_cases', yll_cons, f.values)]:
    ex = float(((s.d1 - s.E)*w).sum()); yll = float(((s.d1 - s.E)*w*s.ypd).sum())
    exs = (ex_sim*w).sum(1); daly = gbd40['yld'] + yll; dsim = gbd40['yld'] + sim
    res[lab] = dict(excess=ex, excess_ci=[float(np.quantile(exs, .025)), float(np.quantile(exs, .975))], ratio_vs_gbd=ex/gbd40['deaths'],
                    excess_yll=yll, excess_yll_ci=[float(np.quantile(sim, .025)), float(np.quantile(sim, .975))],
                    daly=daly, daly_ratio=daly/gbd40['daly'], daly_ratio_ci=[float(np.quantile(dsim/gbd40['daly'], .025)), float(np.quantile(dsim/gbd40['daly'], .975))],
                    yll_share=yll/daly)
res['gbd_ge40']['yll_share'] = gbd40['yll']/gbd40['daly']
json.dump(res, open(out, 'w'), indent=1); print(json.dumps(res, indent=1))
