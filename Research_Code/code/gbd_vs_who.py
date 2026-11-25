"""GBD 2023 RA deaths (ages >=20) versus WHO-registered underlying-cause M05+M06 deaths, 2019 and 2023; plus registered
counts used as the certified comparator for France 2019 and Denmark 2000-18. Args: who_csv gbd_dir t5_csv out_csv"""
import sys, io, zipfile
import numpy as np, pandas as pd
who, gbddir, t5p, out = sys.argv[1:5]
w = pd.read_csv(who)
def gz(name):
    z = zipfile.ZipFile(f'{gbddir}/{name}'); m = [m for m in z.infolist() if m.filename.endswith('.csv')][0]
    return pd.read_csv(io.BytesIO(z.read(m)), low_memory=False)
g = pd.concat([gz('IHME-GBD_2023_DATA-aaf6002d-1.zip'), gz('IHME-GBD_2023_DATA-c82898aa-1.zip')])
g = g[(g.metric_name == 'Number') & (g.measure_name == 'Deaths') & g.age_name.eq('20+ years') & g.sex_name.eq('Both')]
rows = []
for yr in [2019, 2023]:
    for _, r in w[w.Year == yr].iterrows():
        q = g[(g.location_name == r.country) & (g.year == yr)]
        if len(q): rows.append(dict(country=r.country, year=yr, registered=r.M05_M06, gbd=float(q.val.iloc[0]), gbd_lo=float(q.lower.iloc[0]), gbd_hi=float(q.upper.iloc[0])))
c = pd.DataFrame(rows); c['gbd_to_registered'] = c.gbd/c.registered; c.to_csv(out, index=False)
print(c.round(2).to_string(index=False))
print('median GBD/registered 2019: %.3f; range %.3f-%.3f' % (c[c.year == 2019].gbd_to_registered.median(), c[c.year == 2019].gbd_to_registered.min(), c[c.year == 2019].gbd_to_registered.max()))
dk = w[(w.country == 'Denmark') & w.Year.between(2000, 2018)]; fr = w[(w.country == 'France') & (w.Year == 2019)]
t5 = pd.read_csv(t5p)
print('Denmark registered M05+M06 2000-18 mean %.1f per year (M05 %.1f, M06 %.1f); n years %d' % (dk.M05_M06.mean(), dk.M05.mean(), dk.M06.mean(), len(dk)))
print('France registered 2019 M05+M06 %d' % fr.M05_M06.iloc[0])
rng = np.random.default_rng(20261123); B = 20000
for _, r in t5.iterrows():
    if 'France' in r.setting: C = float(fr.M05_M06.iloc[0]); years = 1
    elif 'Denmark' in r.setting: C = float(dk.M05_M06.sum()); years = len(dk)
    else: continue
    ex = rng.normal(r.excess, (r.excess_hi - r.excess_lo)/(2*1.959964), B); cc = np.maximum(rng.poisson(C, B), 1)/years
    print('%s: excess %.0f per year / registered %.1f per year = %.2f (%.2f-%.2f)' % (r.setting, r.excess, C/years, r.excess/(C/years), *np.quantile(ex/cc, [.025, .975])))
