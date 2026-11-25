"""Main analysis: excess all-cause deaths among people with rheumatoid arthritis (RA) versus deaths
attributed to RA. Inputs are read as data; paths are passed as arguments.
  1 ccdss   CCDSS 2025 RA crude all-cause mortality table (Canada excl. New Brunswick; 5-year groups)
  2 panel   provincial panel (prepare/build_provincial_panel.py; ages 50-64, 65-79, >=80 by sex)
  3 statcan Statistics Canada table 13-10-0150-01 zip (underlying cause of death, ICD-10 M05-M06)
  4 gbddir  GBD 2023 results-tool files for 2019 and 2023 (RA and all-cause measures, 204 countries and the world)
  5 cells   GBD cells with the four mortality profiles (prepare/build_gbd_cells.py; ages 45-89)
  6 frdir   French SNDS table 2 + Eurostat 2019 deaths and populations
  7 atlas   Danish Atlas text files
  8 out     output folder
Excess deaths = observed deaths among people with RA minus deaths expected at the rates of the whole
population (people with and without RA) in the same sex and 5-year age group and year."""
import sys, io, json, zipfile
from pathlib import Path
import numpy as np, pandas as pd
ccdss, panel, statcan, gbddir, cells_path, frdir, atlasdir, out = sys.argv[1:9]
out = Path(out); out.mkdir(parents=True, exist_ok=True)
Z = 1.959964; rng = np.random.default_rng(20261123); B = 10000; K = {}
# ------------------------------------------------------------------ CCDSS national
raw = open(ccdss, encoding='latin-1').read().splitlines()
t = pd.read_csv(io.StringIO('\n'.join(raw[1:])), low_memory=False)
t.columns = [c.replace('\xa0', ' ') for c in t.columns]
for c in ['Age group', 'Fiscal year']: t[c] = t[c].astype(str).str.replace('\x96', '–').str.replace('*', '', regex=False)
nm = lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')
t['d1'] = nm(t['Counts of mortality with the condition']); t['n1'] = nm(t['Population of mortality with the condition'])
t['d0'] = nm(t['Counts of mortality without the condition']); t['n0'] = nm(t['Population of mortality without the condition'])
c5 = t[(t.Geography == 'Canada') & (t['Age group type'] == '5 Year') & t.Sex.isin(['Females', 'Males'])].copy()
c5['age'] = c5['Age group'].str.extract(r'^(\d+)').astype(int); c5['year'] = c5['Fiscal year'].str[:4].astype(int)
c5['sex'] = c5.Sex.map({'Females': 'Women', 'Males': 'Men'})
c5 = c5[c5.age >= 40].copy()
assert c5[['d1', 'n1', 'd0', 'n0']].notna().all().all(), 'cells at ages >=40 must be complete'
assert c5.groupby('year').size().eq(22).all()          # 11 age groups (40-44 ... 90+) x 2 sexes
c5['N'] = c5.n1 + c5.n0; c5['mpop'] = (c5.d1 + c5.d0) / c5.N
c5['E'] = c5.n1 * c5.mpop; c5['E0'] = c5.n1 * c5.d0 / c5.n0
w0 = c5[c5.year.between(2007, 2009)].groupby(['age', 'sex']).n1.sum(); w0 = (w0 / w0.sum()).rename('w')   # 2007-09 RA case standard
c5 = c5.merge(w0.reset_index(), on=['age', 'sex'])
# ------------------------------------------------------------------ Statistics Canada coded deaths
sc = pd.read_csv(io.BytesIO(zipfile.ZipFile(statcan).read('13100150.csv')), low_memory=False)
sc = sc[sc['Cause of death (ICD-10)'].isin(['Seropositive rheumatoid arthritis [M05]', 'Other rheumatoid arthritis [M06]'])]
AGEMAP = {f'{a} to {a+4} years': a for a in range(0, 90, 5)} | {'90 years and over': 90}
sc['age'] = sc['Age group'].map(AGEMAP); sc['sex'] = sc.Sex.map({'Females': 'Women', 'Males': 'Men', 'Both sexes': 'Both'})
coded = sc[sc.age.notna()].groupby(['REF_DATE', 'sex', 'age']).VALUE.sum()
def coded_n(year, lo=40, hi=200, sex='Both'):
    return float(sum(coded.get((year, sex, a), 0) for a in list(range(0, 90, 5)) + [90] if lo <= a < hi))
coded_all_ages = sc[sc['Age group'].eq('Total, all ages') & sc.sex.eq('Both')].groupby('REF_DATE').VALUE.sum()
# ------------------------------------------------------------------ GBD exports
FIVE = {f'{a}-{a+4} years': a for a in range(20, 95, 5)} | {'95+ years': 95}
def gz(name):
    z = zipfile.ZipFile(f'{gbddir}/{name}'); m = [m for m in z.infolist() if m.filename.endswith('.csv')][0]
    return pd.read_csv(io.BytesIO(z.read(m)), low_memory=False)
gra = pd.concat([gz('IHME-GBD_2023_DATA-aaf6002d-1.zip'), gz('IHME-GBD_2023_DATA-c82898aa-1.zip')])
gnum = gra[gra.metric_name.eq('Number')]
def gbd(loc, yr, measure, lo=None, hi=None, agg=None, col='val'):
    d = gnum[(gnum.location_name == loc) & (gnum.year == yr) & (gnum.measure_name == measure)]
    if agg: return float(d[(d.age_name == agg) & d.sex_name.eq('Both')][col].sum())
    d = d[d.sex_name.isin(['Female', 'Male']) & d.age_name.isin(FIVE)].copy(); d['age'] = d.age_name.map(FIVE)
    return float(d[(d.age >= lo) & (d.age < hi)][col].sum())
# ------------------------------------------------------------------ helpers
def excess_row(s):
    D = s.d1.sum(); E = s.E.sum(); ex = D - E
    return dict(people_ra=s.n1.sum(), deaths=D, expected=E, excess=ex, excess_lo=ex - Z*np.sqrt(D), excess_hi=ex + Z*np.sqrt(D),
                smr=D/E, smr_lo=(D - Z*np.sqrt(D))/E, smr_hi=(D + Z*np.sqrt(D))/E,
                excess_vs_nonra=D - s.E0.sum(), excess_per_1000=1000*ex/s.n1.sum(), prevalence_pct=100*s.n1.sum()/s.N.sum())
def ratio_ci(D, E, C):
    sim = (rng.poisson(D, B) - E) / np.maximum(rng.poisson(C, B), 1)
    return float(np.quantile(sim, 0.025)), float(np.quantile(sim, 0.975))
def std_rates(g):
    w = g.w / g.w.sum(); m1 = g.d1/g.n1
    pr = (w*m1).sum(); po = (w*g.mpop).sum()
    vpr = (w**2*g.d1/g.n1**2).sum(); vpo = (w**2*(g.d1+g.d0)/g.N**2).sum(); cov = (w**2*g.d1/(g.n1*g.N)).sum()
    l = np.log(pr/po); se = np.sqrt(vpr/pr**2 + vpo/po**2 - 2*cov/(pr*po))
    return dict(patient_rate=1000*pr, population_rate=1000*po, gap=1000*(pr-po), srr=np.exp(l), srr_lo=np.exp(l-Z*se), srr_hi=np.exp(l+Z*se), log_srr=l, se=se)
AGEB = {'40–64': (40, 65), '65–79': (65, 80), '≥80': (80, 200), 'All ≥40': (40, 200)}
# ------------------------------------------------------------------ 1. annual series, ages >=40
rows = []
for yr, s in c5.groupby('year'):
    r = excess_row(s) | dict(year=yr, coded_statcan=coded_n(yr), coded_statcan_all_ages=float(coded_all_ages.get(yr, np.nan)))
    r['ratio'] = r['excess'] / r['coded_statcan']; r['ratio_lo'], r['ratio_hi'] = ratio_ci(r['deaths'], r['expected'], r['coded_statcan'])
    r['coded_pct_of_excess'] = 100*r['coded_statcan']/r['excess']; r['coded_pct_of_deaths'] = 100*r['coded_statcan']/r['deaths']
    rows.append(r)
annual = pd.DataFrame(rows); annual.to_csv(out/'T1_canada_annual_excess_vs_coded.csv', index=False)
# ------------------------------------------------------------------ 2. by age and sex, 2019 and 2023
rows = []
for yr in [2019, 2023]:
    for sex in ['Both', 'Women', 'Men']:
        for b, (lo, hi) in AGEB.items():
            s = c5[(c5.year == yr) & (c5.age >= lo) & (c5.age < hi) & ((c5.sex == sex) if sex != 'Both' else True)]
            r = excess_row(s) | dict(year=yr, sex=sex, ages=b, coded_statcan=coded_n(yr, lo, hi, sex))
            if sex == 'Both':
                r['gbd_coded'] = gbd('Canada', yr, 'Deaths', lo, hi); r['gbd_people_ra'] = gbd('Canada', yr, 'Prevalence', lo, hi)
            r['ratio'] = r['excess']/r['coded_statcan'] if r['coded_statcan'] > 0 else np.nan
            if r['coded_statcan'] > 0: r['ratio_lo'], r['ratio_hi'] = ratio_ci(r['deaths'], r['expected'], r['coded_statcan'])
            rows.append(r)
byage = pd.DataFrame(rows)
byage['conservative_excess'] = byage.excess / byage.people_ra * byage.gbd_people_ra      # excess per CCDSS case x GBD case count
byage['conservative_ratio'] = byage.conservative_excess / byage.coded_statcan
byage.to_csv(out/'T2_canada_by_age_sex_2019_2023.csv', index=False)
# ------------------------------------------------------------------ 3. relative excess over time (direct standardisation)
rows = []
for b, (lo, hi) in AGEB.items():
    g = c5[(c5.age >= lo) & (c5.age < hi)]
    for yr, s in g.groupby('year'): rows.append(std_rates(s) | dict(ages=b, year=yr))
srr = pd.DataFrame(rows); srr.to_csv(out/'T3_canada_standardised_annual.csv', index=False)
tr = []
for b in AGEB:
    a = srr[(srr.ages == b) & srr.year.between(2010, 2019)]
    for col, se_col, scale in [('log_srr', 'se', 100)]:
        W = 1/a[se_col].values**2; X = np.c_[np.ones(len(a)), a.year - 2010]; A = X.T@(X*W[:, None])
        beta = np.linalg.solve(A, X.T@(W*a[col].values)); res = a[col].values - X@beta
        disp = max(1, (W*res**2).sum()/(len(a)-2)); sb = np.sqrt(np.linalg.inv(A)[1, 1]*disp)
        tr.append(dict(ages=b, measure='SRR % change per year, 2010-19', est=scale*beta[1], lo=scale*(beta[1]-Z*sb), hi=scale*(beta[1]+Z*sb)))
    for col in ['gap', 'patient_rate', 'population_rate']:
        sd = {'gap': None}
        v = a[col].values; X = np.c_[np.ones(len(a)), a.year - 2010]; beta = np.linalg.lstsq(X, v, rcond=None)[0]
        tr.append(dict(ages=b, measure=f'{col} change per year 2010-19 (OLS)', est=beta[1]))
trend = pd.DataFrame(tr); trend.to_csv(out/'T3b_trend_2010_19.csv', index=False)
# pandemic: pooled SRR vs 2015-19 (bootstrap on cell deaths)
def pooled(g, yrs):
    q = g[g.year.isin(yrs)].groupby(['age', 'sex'], as_index=False).agg(d1=('d1', 'sum'), n1=('n1', 'sum'), d0=('d0', 'sum'), n0=('n0', 'sum'), w=('w', 'first'))
    q['N'] = q.n1 + q.n0; q['mpop'] = (q.d1 + q.d0)/q.N; return q
def boot_srr(q, n=4000):
    D1 = rng.poisson(q.d1.values, (n, len(q))); D0 = rng.poisson(q.d0.values, (n, len(q))); w = q.w.values/q.w.sum()
    return (w*D1/q.n1.values).sum(1) / (w*(D1+D0)/q.N.values).sum(1)
rows = []
for b, (lo, hi) in AGEB.items():
    g = c5[(c5.age >= lo) & (c5.age < hi)]; ref = pooled(g, range(2015, 2020)); r0 = std_rates(ref)['srr']; br0 = boot_srr(ref)
    for lab, yrs in [('2020', [2020]), ('2021', [2021]), ('2022', [2022]), ('2023', [2023]), ('2021–22', [2021, 2022])]:
        q = pooled(g, yrs); r1 = std_rates(q)['srr']; ratio = boot_srr(q)/br0
        rows.append(dict(ages=b, period=lab, srr_2015_19=r0, srr=r1, ratio=r1/r0, lo=float(np.exp(np.quantile(np.log(ratio), .025))), hi=float(np.exp(np.quantile(np.log(ratio), .975)))))
pand = pd.DataFrame(rows); pand.to_csv(out/'T3c_pandemic_vs_2015_19.csv', index=False)
# ------------------------------------------------------------------ 4. provinces, ages >=50, pooled complete province-years 2010-23
p = pd.read_csv(panel); p = p[p.geography != 'Canada']
cnt = p.groupby(['geography', 'year']).size(); full = cnt[cnt == 6].index
p = p.set_index(['geography', 'year']).loc[full].reset_index(); p['N'] = p.n_ra + p.n_nonra
p['mpop'] = (p.d_ra + p.d_nonra)/p.N; p['E'] = p.n_ra*p.mpop
popw = p.groupby(['age_group', 'sex']).N.sum(); popw = popw/popw.sum()
rw = p.groupby(['age_group', 'sex']).weight.first(); rw = rw/rw.sum()
prov = []
for geo, gq in p.groupby('geography'):
    c = gq.groupby(['age_group', 'sex'])[['d_ra', 'n_ra', 'd_nonra', 'n_nonra', 'N', 'E']].sum()
    m1 = c.d_ra/c.n_ra; mp = (c.d_ra + c.d_nonra)/c.N
    pr = (rw*m1).sum(); po = (rw*mp).sum(); vpr = (rw**2*c.d_ra/c.n_ra**2).sum(); vpo = (rw**2*(c.d_ra+c.d_nonra)/c.N**2).sum(); cov = (rw**2*c.d_ra/(c.n_ra*c.N)).sum()
    l = np.log(pr/po); sl = np.sqrt(vpr/pr**2 + vpo/po**2 - 2*cov/(pr*po))
    exr = (c.d_ra - c.E)/c.N; pe = (popw*exr).sum()*1e5; spe = np.sqrt((popw**2*c.d_ra/c.N**2).sum())*1e5
    prov.append(dict(province=geo, years=gq.year.nunique(), first=int(gq.year.min()), last=int(gq.year.max()), deaths=c.d_ra.sum(), expected=c.E.sum(),
        excess=(c.d_ra - c.E).sum(), srr=np.exp(l), srr_lo=np.exp(l - Z*sl), srr_hi=np.exp(l + Z*sl),
        excess_per_100k=pe, excess_per_100k_lo=pe - Z*spe, excess_per_100k_hi=pe + Z*spe, prevalence_std_pct=(popw*c.n_ra/c.N).sum()*100,
        mean_people_ra=gq.groupby('year').n_ra.sum().mean()))
prov = pd.DataFrame(prov).sort_values('excess_per_100k'); prov.to_csv(out/'T4_provinces_50plus.csv', index=False)
def het(est, lo, hi, log=False):
    e = np.log(est) if log else est; se = ((np.log(hi) - np.log(lo)) if log else (hi - lo))/(2*Z); w = 1/se**2
    mu = (w*e).sum()/w.sum(); Q = float((w*(e - mu)**2).sum()); return dict(Q=Q, df=len(e)-1, I2=max(0, (Q-(len(e)-1))/Q))
K['provinces'] = dict(srr=het(prov.srr, prov.srr_lo, prov.srr_hi, True), pop=het(prov.excess_per_100k, prov.excess_per_100k_lo, prov.excess_per_100k_hi),
    spearman_prev_srr=float(prov.prevalence_std_pct.rank().corr(prov.srr.rank())), spearman_prev_pop=float(prov.prevalence_std_pct.rank().corr(prov.excess_per_100k.rank())))
early = p[p.year <= 2015].groupby('geography').apply(lambda g: ((g.d_ra - g.E).sum()), include_groups=False)
# ------------------------------------------------------------------ 5. France 2019, Denmark 2000-18, Veneto 2012-20
t2 = pd.read_csv(f'{frdir}/france_snds_2019_table2.csv')
def ej(path):
    j = json.load(open(path)); dims = j['id']; size = j['size']; cat = {k: j['dimension'][k]['category']['index'] for k in dims}
    keys = {k: (sorted(cat[k], key=lambda x: cat[k][x]) if isinstance(cat[k], dict) else cat[k]) for k in dims}; rows = []
    for flat, v in j['value'].items():
        i = int(flat); idx = []
        for s in reversed(size): idx.append(i % s); i //= s
        idx = idx[::-1]; rows.append({k: keys[k][idx[n]] for n, k in enumerate(dims)} | {'value': v})
    return pd.DataFrame(rows)
def single(a):
    if a == 'Y_LT1': return 0
    if a.startswith('Y') and a[1:].isdigit(): return int(a[1:])
    return 100 if a == 'Y_OPEN' else None
fd = ej(f'{frdir}/eurostat_fr_deaths_2019.json'); fp = ej(f'{frdir}/eurostat_fr_population_2019_2020.json')
for d in (fd, fp): d['a'] = d['age'].map(single)
fd = fd[fd.a.notna() & fd.sex.isin(['M', 'F'])]; fp = fp[fp.a.notna() & fp.sex.isin(['M', 'F'])]
fr = fd.groupby(['sex', 'a']).value.sum().rename('deaths').reset_index().merge(fp.groupby(['sex', 'a']).value.mean().rename('pop').reset_index(), on=['sex', 'a'])
fr['age_start'] = (fr.a//5*5).clip(upper=95); fr5 = fr.groupby(['sex', 'age_start'])[['deaths', 'pop']].sum().reset_index(); fr5['rate'] = fr5.deaths/fr5['pop']
long = pd.concat([t2[['age_start', 'female_RA_patients']].rename(columns={'female_RA_patients': 'n'}).assign(sex='F'),
                  t2[['age_start', 'male_RA_patients']].rename(columns={'male_RA_patients': 'n'}).assign(sex='M')]).merge(fr5, on=['sex', 'age_start'])
E_fr = float((long.n*long.rate).sum()); D_fr = 9722.0
intl = [dict(setting='France (SNDS), 2019', case_definition='Hospital diagnosis or long-term illness status, M05/M06', deaths=D_fr, expected=E_fr,
             excess=D_fr - E_fr, excess_lo=D_fr - E_fr - Z*np.sqrt(D_fr), excess_hi=D_fr - E_fr + Z*np.sqrt(D_fr), smr=D_fr/E_fr,
             coded=gbd('France', 2019, 'Deaths', agg='20+ years'), coded_lo=gbd('France', 2019, 'Deaths', agg='20+ years', col='lower'),
             coded_hi=gbd('France', 2019, 'Deaths', agg='20+ years', col='upper'), coded_source='GBD 2023, 2019, ages ≥20', people=float(long.n.sum()))]
me = pd.read_csv(f'{atlasdir}/atlas_02_main_estimates.txt', sep=' ', quotechar='"'); mr = pd.read_csv(f'{atlasdir}/atlas_04_MRR.txt', sep=' ', quotechar='"')
dk_c = gbd('Denmark', 2019, 'Deaths', agg='20+ years')
for code, lab in [(10919, 'M05'), (10920, 'M06')]:
    D = float(me[(me.id == code) & (me.sex == 'All')].n_deaths.iloc[0]); r = mr[(mr.id == code) & (mr.sex == 'All') & (mr.model == 'All')].iloc[0]
    intl.append(dict(setting=f'Denmark (Atlas) {lab}, 2000–18, annual mean', case_definition=f'Hospital contact with {lab}', deaths=D/19, expected=D/19/r.MRR,
        excess=D*(1 - 1/r.MRR)/19, excess_lo=D*(1 - 1/r.MRR_left)/19, excess_hi=D*(1 - 1/r.MRR_right)/19, smr=r.MRR, coded=dk_c,
        coded_lo=gbd('Denmark', 2019, 'Deaths', agg='20+ years', col='lower'), coded_hi=gbd('Denmark', 2019, 'Deaths', agg='20+ years', col='upper'),
        coded_source='GBD 2023, 2019, ages ≥20', people=float(me[(me.id == code) & (me.sex == 'All')].prev_cases.iloc[0])))
# Veneto: same death records; SMR 1.28 (1.21-1.35) from incident cohort; RA underlying cause in 355 of 7435 deaths
V_D, V_C = 7435.0, 355.0
intl.append(dict(setting='Veneto, Italy, 2012–20 (prevalent cohort)', case_definition='Co-payment exemption or hospital diagnosis', deaths=V_D,
    expected=V_D/1.28, excess=V_D*(1 - 1/1.28), excess_lo=V_D*(1 - 1/1.21), excess_hi=V_D*(1 - 1/1.35), smr=1.28, coded=V_C, coded_lo=np.nan, coded_hi=np.nan,
    coded_source='Underlying cause on the same death records', people=37996.0))
intl = pd.DataFrame(intl); intl['ratio'] = intl.excess/intl.coded
intl['ratio_lo'] = intl.excess_lo/intl.coded; intl['ratio_hi'] = intl.excess_hi/intl.coded
intl.loc[intl.coded_lo.notna(), 'ratio_lo_inclUI'] = intl.excess_lo/intl.coded_hi; intl.loc[intl.coded_lo.notna(), 'ratio_hi_inclUI'] = intl.excess_hi/intl.coded_lo
intl.to_csv(out/'T5_international.csv', index=False)
# ------------------------------------------------------------------ 6. global scenarios, ages 45-89
cl = pd.read_csv(cells_path, low_memory=False)
cl = cl[(cl.location_name != 'Global') & (cl.scenario != 'Canada 2019 fixed')].copy(); cl['age'] = cl.age_name.str.extract(r'^(\d+)').astype(int)
ypd = (cl.groupby(['sex_name', 'age', 'year']).ylls.sum()/cl.groupby(['sex_name', 'age', 'year']).ra_deaths.sum()).rename('ypd').reset_index()
cl = cl.merge(ypd, on=['sex_name', 'age', 'year']); cl['excess_yll'] = cl.equivalent*cl.ypd
base = cl[cl.scenario == 'Canada relative']; rows = []
for yr in [2019, 2023]:
    b = base[base.year == yr]
    for s_, d in cl[cl.year == yr].groupby('scenario'):
        rows.append(dict(year=yr, scenario=s_, cases=b.cases.sum(), gbd_coded=b.ra_deaths.sum(), modelled_excess=d.equivalent.sum(),
            ratio=d.equivalent.sum()/b.ra_deaths.sum(), gbd_yll=b.ylls.sum(), gbd_yld=b.ylds.sum(), excess_yll=d.excess_yll.sum(),
            daly_gbd=b.ylls.sum() + b.ylds.sum(), daly_with_excess=b.ylds.sum() + d.excess_yll.sum(),
            share_80_89=d[d.age >= 80].equivalent.sum()/d.equivalent.sum()))
glob = pd.DataFrame(rows); glob['daly_ratio'] = glob.daly_with_excess/glob.daly_gbd; glob['yll_share_with_excess'] = glob.excess_yll/glob.daly_with_excess
glob.to_csv(out/'T6_global_45_89.csv', index=False)
can_scen = cl[cl.location_name == 'Canada'].groupby(['scenario', 'year']).agg(model_excess=('equivalent', 'sum'), gbd_cases=('cases', 'sum'), gbd_coded=('ra_deaths', 'sum')).reset_index()
o = c5[(c5.age >= 45) & (c5.age < 90)].groupby('year').apply(lambda s: pd.Series(dict(obs_excess=s.d1.sum() - s.E.sum(), ccdss_people=s.n1.sum())), include_groups=False).reset_index()
can_scen = can_scen.merge(o, on='year'); can_scen.to_csv(out/'T6b_canada_scenario_vs_observed_45_89.csv', index=False)
# ------------------------------------------------------------------ key numbers
def f(x): return float(x) if isinstance(x, (np.floating, np.integer, int, float)) else x
K['annual_2023'] = {k: f(v) for k, v in annual[annual.year == 2023].iloc[0].items()}
K['annual_2019'] = {k: f(v) for k, v in annual[annual.year == 2019].iloc[0].items()}
K['annual_range'] = dict(ratio_min=f(annual.ratio.min()), ratio_max=f(annual.ratio.max()), coded_pct_excess_min=f(annual.coded_pct_of_excess.min()),
    coded_pct_excess_max=f(annual.coded_pct_of_excess.max()), years=f'{annual.year.min()}-{annual.year.max()}', total_excess=f(annual.excess.sum()), total_coded=f(annual.coded_statcan.sum()))
K['gbd_canada_2023_ge40'] = gbd('Canada', 2023, 'Deaths', 40, 200); K['gbd_canada_2019_ge40'] = gbd('Canada', 2019, 'Deaths', 40, 200)
K['gbd_canada_2023_20plus'] = [gbd('Canada', 2023, 'Deaths', agg='20+ years', col=c) for c in ['val', 'lower', 'upper']]
K['gbd_canada_people_2023_ge40'] = gbd('Canada', 2023, 'Prevalence', 40, 200)
K['gbd_global_20plus'] = {str(y): [gbd('Global', y, 'Deaths', agg='20+ years', col=c) for c in ['val', 'lower', 'upper']] for y in [2019, 2023]}
json.dump(K, open(out/'key_numbers.json', 'w'), indent=1, ensure_ascii=False, default=f)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
print(annual[['year', 'people_ra', 'deaths', 'expected', 'excess', 'excess_lo', 'excess_hi', 'smr', 'excess_per_1000', 'coded_statcan', 'ratio', 'ratio_lo', 'ratio_hi', 'coded_pct_of_excess']].round(2).to_string(index=False))
print(byage[['year', 'sex', 'ages', 'people_ra', 'deaths', 'expected', 'excess', 'excess_lo', 'excess_hi', 'smr', 'excess_per_1000', 'coded_statcan', 'gbd_coded', 'ratio', 'ratio_lo', 'ratio_hi', 'conservative_excess', 'conservative_ratio']].round(2).to_string(index=False))
print(srr[srr.year.isin([2010, 2015, 2019, 2020, 2021, 2022, 2023])][['ages', 'year', 'patient_rate', 'population_rate', 'gap', 'srr', 'srr_lo', 'srr_hi']].round(3).to_string(index=False))
print(trend.round(3).to_string(index=False)); print(pand.round(3).to_string(index=False)); print(prov.round(2).to_string(index=False))
print(intl[['setting', 'deaths', 'expected', 'excess', 'excess_lo', 'excess_hi', 'smr', 'coded', 'coded_lo', 'coded_hi', 'ratio', 'ratio_lo', 'ratio_hi', 'ratio_lo_inclUI', 'ratio_hi_inclUI']].round(2).to_string(index=False))
print(glob.round(3).to_string(index=False)); print(can_scen.round(1).to_string(index=False))
print(json.dumps({k: K[k] for k in ['provinces', 'annual_range', 'gbd_canada_2023_ge40', 'gbd_canada_2019_ge40', 'gbd_canada_2023_20plus', 'gbd_canada_people_2023_ge40', 'gbd_global_20plus']}, indent=1, default=f))
