# Analysis code: excess mortality in people with rheumatoid arthritis and its under-representation in cause-of-death statistics

This repository contains the statistical analysis code for "Excess mortality in people with rheumatoid arthritis and its under-representation in cause-of-death statistics: a multinational population-based study" (Zhao Z, et al.). The code covers data preparation, variable construction, all estimates, and the consistency checks used in the paper. It produces every number, table value, and plotted series from the raw inputs.

The data are distributed separately as Research_Data, the supplementary data file submitted with the manuscript. That package holds the open raw inputs, the analysis datasets, the published results, a data dictionary, and download instructions for third-party files that cannot be redistributed.

## What the analysis does

1. **Canada (CCDSS, ages ≥40 years, fiscal years 2007–08 to 2023–24).**
   - Expected deaths among people with RA = people with RA × whole-population mortality in the same year, sex, and 5-year age group.
   - Excess deaths = observed − expected.
   - Excess deaths are compared with deaths with ICD-10 M05–M06 as the underlying cause (Statistics Canada).
   - 95% CIs come from Poisson simulation with a fixed seed (20261123).
   - Rates are directly standardised to the 2007–09 age–sex distribution of people with RA. Changes are compared between 2007–09, 2017–19, and 2021–23.
   - Conservative analysis: excess deaths per person are applied to GBD 2023 case numbers.
   - Further analyses: provinces, the pandemic period, and a comparator of people without RA.
2. **Europe.** The same comparison is made for France (SNDS, 2019), Denmark (national mortality atlas, M05 and M06, 2000–18), and Veneto, Italy (2012–20), using published counts and registered deaths from the WHO Mortality Database.
3. **GBD 2023 trends (1990–2023).**
   - Observed GBD RA deaths per person with RA are compared with the values expected had they changed in proportion to all-cause mortality.
   - Locations are weighted by their RA cases.
   - Rates are standardised to the global 2019 RA age–sex distribution.
   - The comparison is made for the world, the seven super-regions, 11 high-income countries, and 204 countries.
   - GBD RA deaths are also compared with registered deaths.
4. **Global excess deaths and burden (ages 45–89 years, 2019 and 2023).**
   - Four mortality profiles are applied to GBD 2023 populations: Canadian relative and additive profiles, and Danish M05 and M06 profiles.
   - Excess years of life lost are added to GBD years lived with disability, giving recomposed DALYs.


## Requirements

- Python 3.14 (tested; Python 3.11 or later should work).
- pandas 3.0.1 and NumPy 2.4.2 (`pip install -r requirements.txt`).
- zsh, for `run_all.sh` (default shell on macOS; available on Linux).
- No other software and no internet access are needed at run time.

Every script is run with `python3 -I` (isolated mode). Input paths are passed as arguments. The scripts read input files only as data.

## Data layout

Place the data package next to this repository (clone or unpack Research_Code, and unpack Research_Data beside it):

```
Research_Code/        this repository
Research_Data/        data package
  raw_open/           CCDSS, Statistics Canada, Eurostat, and French SNDS inputs (included)
  raw_restricted/     GBD 2023, WHO Mortality Database, and Danish atlas files (download; see DOWNLOAD_INSTRUCTIONS.md)
  analysis_data/      provincial panel, GBD location map, WHO extract
  results/            published outputs used in the manuscript
```

All paths are set in `config.sh`, relative to this folder. To use another location, set `DATA`, `RESTRICTED`, or a single input variable (for example `GBD_ANNUAL=/path/to/files ./run_all.sh`).

## Running the analysis

```
./run_all.sh output
python3 -I check_against_published.py output ../Research_Data/results
```

`run_all.sh` runs the steps below in order. It stops at the first error and writes one log per step to `output/logs/`. On a laptop it takes about 70 seconds.

| Step | Script | Inputs | Outputs |
|---|---|---|---|
| 1 | `prepare/build_provincial_panel.py` | CCDSS export | `prepared/provincial_panel.csv` |
| 2 | `prepare/build_gbd_cells.py` | GBD 2019/2023 files, CCDSS export, Danish atlas age-specific file | `prepared/gbd_cells_profiles.csv.gz` (GBD cells with the four profiles) |
| 3 | `prepare/build_location_map.py` | GBD tool hierarchy and labels (or the supplied map) | `prepared/gbd_location_map.csv` |
| 4 | `who_ra_deaths.py` | WHO raw files (or the supplied extract) | `prepared/who_ra_underlying_deaths.csv` |
| 5 | `final_analysis.py` | CCDSS, panel, Statistics Canada, GBD 2019/2023, cells, France and Eurostat, Danish atlas | `final/T1`–`T6b`, `key_numbers.json` |
| 6 | `gbd_superregion.py` | location map, cells | `final/T7_superregion_2023.csv` |
| 7 | `std_certified.py`, `std_change_periods.py` | CCDSS, Statistics Canada | standardised series and period changes (`display/`) |
| 8 | `conservative_5yr.py`, `conservative_ci.py` | CCDSS, Statistics Canada, GBD | conservative estimate and its CI |
| 9 | `canada_yll_trend.py`, `canada_daly_table.py` | CCDSS, Statistics Canada, GBD | Canadian YLLs and DALYs |
| 10 | `figure_data.py`, `table_values.py`, `display_tables.py`, `last_numbers.py`, `extra_numbers.py`, `age_shift.py` | outputs above | numbers for the text, table, and appendix tables |
| 11 | `gbd_vs_who.py` | WHO extract, GBD | table S9 |
| 12 | `gbd_trends.py`, `gbd_annual.py` | GBD annual and 2019/2023 files, location map, WHO extract | `gbd_trends/G1`–`G7` |
| 13 | `gbd_expected.py`, `gbd_counterfactual.py`, `gbd_expected_sensitivity.py` | GBD annual files, location map | observed vs expected trends (`G10`–`G12`; figure 3; table S11) |
| 14 | `gbd_trend_checks.py` | GBD files, Statistics Canada | consistency checks of the GBD files (log only) |
| 15 | `gbd_trend_figdata.py`, `figure1b_population.py`, `derive_figure_data.py` | outputs above | numeric series for figures 1–4 and S1 |

The data package README maps each output file to the tables, figures, and text it supports, and `Research_Data/Data_dictionary.csv` defines every field.

## Reproduction check

The code was run on a clean copy of the two packages, with the restricted files placed as described. Results:
- 48 of 49 output files were byte-for-byte identical to the published files in `Research_Data/results`.
- The remaining file (`figdata/Figure1B_standardised_rates.csv`) agreed to within 3 × 10⁻¹⁴ (floating-point rounding).
- All 146 numbers checked in the manuscript, appendix, and covering letter were identical.
- The three analysis datasets supplied in `Research_Data/analysis_data` were rebuilt byte for byte (`output/prepared`).

If a provider has revised its data since the study downloads, re-downloaded files can give slightly different results; the checksums in `DOWNLOAD_INSTRUCTIONS.md` identify the files used.

## Licence

**[TO CONFIRM by the authors and their institution]** Proposed licence: MIT for the code in this repository, which was written for this study. The repository contains no third-party code. Data are not covered by the code licence; see `Research_Data/README.md`.

