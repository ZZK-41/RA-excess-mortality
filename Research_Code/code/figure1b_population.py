"""Whole-population mortality at ages >=40, standardised to the 2007-09 RA age-sex distribution (same standard as the excess and attributed
series in Figure 1B), indexed to its 2007-09 mean. Args: T3_csv out_csv"""
import sys
import pandas as pd
t3, out = sys.argv[1:3]
t = pd.read_csv(t3); t = t[t.ages == 'All ≥40'][['year', 'population_rate', 'patient_rate', 'gap']].sort_values('year')
base = t[t.year.between(2007, 2009)].population_rate.mean()
t['population_index'] = 100*t.population_rate/base
t.to_csv(out, index=False)
chg = 100*(t[t.year.between(2021, 2023)].population_rate.mean()/base - 1)
print(t.round(3).to_string(index=False)); print('change 2007-09 to 2021-23: %.2f%%' % chg)
