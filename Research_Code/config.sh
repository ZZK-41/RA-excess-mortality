# Input locations, relative to the folder that contains run_all.sh. Any of them can be overridden as an environment variable.
# Default layout: the Research_Data package unpacked next to Research_Code (../Research_Data). Raw files that cannot be
# redistributed (GBD results-tool files, WHO Mortality Database files, Danish mortality atlas files) are downloaded by the user into
# Research_Data/raw_restricted/ as described in Research_Data/README.md.
DATA=${DATA:-../Research_Data}
RESTRICTED=${RESTRICTED:-$DATA/raw_restricted}
CCDSS=${CCDSS:-$DATA/raw_open/ccdss/ccdss_RA_age_crude_mortality.csv}
STATCAN=${STATCAN:-$DATA/raw_open/statcan/13100150-eng.zip}
FR=${FR:-$DATA/raw_open/france_eurostat}
GBD_HIER=${GBD_HIER:-$RESTRICTED/gbd_tool_metadata/gbd_results_hierarchy.json}
GBD_META=${GBD_META:-$RESTRICTED/gbd_tool_metadata/gbd_results_metadata_en.json}
GBD_2019_2023=${GBD_2019_2023:-$RESTRICTED/gbd2023/files_2019_2023}
GBD_ANNUAL=${GBD_ANNUAL:-$RESTRICTED/gbd2023/files_annual_1990_2023}
WHO_ZIPS=($RESTRICTED/who_mortality_database/morticd10_part{1..6}.zip)
ATLAS=${ATLAS:-$RESTRICTED/danish_atlas}
# Analysis datasets supplied with Research_Data, used when the corresponding raw files are not available:
LOCMAP_SUPPLIED=$DATA/analysis_data/gbd_location_map.csv
WHO_SUPPLIED=$DATA/analysis_data/who_ra_underlying_deaths.csv
