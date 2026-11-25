"""GBD 2023 RA deaths in Canada at ages >=40 years by year (the GBD line in figure 1A).
Args: trends_dir figdata_dir"""
import sys
from pathlib import Path
import pandas as pd
td, fd = Path(sys.argv[1]), Path(sys.argv[2]); fd.mkdir(parents=True, exist_ok=True)
g2 = pd.read_csv(td/'G2_canada_gbd_ge40_annual.csv'); g2.to_csv(fd/'Figure1_gbd_canada_annual.csv', index=False)
print(g2.tail(3).to_string(index=False))
