"""Figure source data derived from analysis outputs (no new estimates).
Args: results_dir out_dir
  Figure1B_standardised_rates.csv  standardised attributed and excess deaths per 1000 people with RA, and indices (2007-09 mean = 100)
  Figure1B_gbd_canada_index.csv    GBD RA deaths per GBD case, Canada, standardised to surveillance 2007-09 weights (from G7)
  Figure4_dalys_world.csv          2023 rows of T6_global_45_89.csv
  FigureS1_superregions.csv        T7_superregion_2023.csv
  TableS4_provinces.csv            T4_provinces_50plus.csv (appendix table S6)"""
import sys
from pathlib import Path
import pandas as pd
r, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
s = pd.read_csv(r/'display/std_certified_vs_excess.csv').set_index('year')
bc, be = s.loc[2007:2009, 'std_certified_per_1000'].mean(), s.loc[2007:2009, 'std_excess_per_1000'].mean()
f = pd.DataFrame({'year': s.index, 'certified_rate': s.std_certified_per_1000.values, 'certified_lo': s.cert_lo.values, 'certified_hi': s.cert_hi.values,
                  'excess_rate': s.std_excess_per_1000.values, 'excess_lo': s.gap_lo.values, 'excess_hi': s.gap_hi.values})
for c, b in [('certified', bc), ('excess', be)]:
    f[f'{c}_index'] = 100*f[f'{c}_rate']/b; f[f'{c}_lo_index'] = 100*f[f'{c}_lo']/b; f[f'{c}_hi_index'] = 100*f[f'{c}_hi']/b
f.to_csv(out/'Figure1B_standardised_rates.csv', index=False)
pd.read_csv(r/'gbd_trends/G7_canada_gbd_ccdss_weights.csv').to_csv(out/'Figure1B_gbd_canada_index.csv', index=False)
t6 = pd.read_csv(r/'final/T6_global_45_89.csv')
t6[t6.year == 2023][['scenario', 'gbd_yld', 'gbd_yll', 'excess_yll', 'daly_gbd', 'daly_with_excess', 'daly_ratio', 'yll_share_with_excess']].to_csv(out/'Figure4_dalys_world.csv', index=False)
pd.read_csv(r/'final/T7_superregion_2023.csv').to_csv(out/'FigureS1_superregions.csv', index=False)
pd.read_csv(r/'final/T4_provinces_50plus.csv').to_csv(out/'TableS4_provinces.csv', index=False)
print('written', out)
