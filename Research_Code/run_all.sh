#!/bin/zsh
# Reproduce every analysis output used in the manuscript and appendix.
# Usage:  ./run_all.sh [output_folder]        (default: ./output)
# Input paths are set in config.sh (relative to this folder); override any of them as environment variables.
set -e
HERE=${0:a:h}; cd $HERE
source ./config.sh
OUT=${1:-output}; S=code; P=code/prepare
for v in CCDSS STATCAN GBD_2019_2023 GBD_ANNUAL FR ATLAS; do [[ -e ${(P)v} ]] || { echo "missing input $v: ${(P)v} (see README and Research_Data/raw_restricted/DOWNLOAD_INSTRUCTIONS.md)"; exit 1; }; done
HAVE_WHO=1; for z in $WHO_ZIPS; do [[ -e $z ]] || HAVE_WHO=0; done
mkdir -p $OUT/prepared $OUT/final $OUT/display $OUT/figdata $OUT/gbd_trends $OUT/logs
run() { local d=$1 name=$2; shift 2; echo "== $name"; python3 -I $d/$name.py "$@" > $OUT/logs/$name.log 2>&1 || { echo "failed: $name (see $OUT/logs/$name.log)"; exit 1; } }
# 1. analysis inputs derived from the raw files
PANEL=$OUT/prepared/provincial_panel.csv; CELLS=$OUT/prepared/gbd_cells_profiles.csv.gz; LOCMAP=$OUT/prepared/gbd_location_map.csv
WHO=$OUT/prepared/who_ra_underlying_deaths.csv
run $P build_provincial_panel   $CCDSS $PANEL
run $P build_gbd_cells          $GBD_2019_2023 $CCDSS $ATLAS/atlas_06_incidences_by_age.txt $CELLS $OUT/prepared/gbd_analysis_cells.csv
if [[ -e $GBD_HIER && -e $GBD_META ]]; then run $P build_location_map $GBD_HIER $GBD_META $GBD_2019_2023 $LOCMAP
else echo "== build_location_map: GBD tool metadata not found; using $LOCMAP_SUPPLIED"; cp $LOCMAP_SUPPLIED $LOCMAP; fi
if (( HAVE_WHO )); then run $S who_ra_deaths $WHO $WHO_ZIPS
else echo "== who_ra_deaths: WHO Mortality Database files not found; using $WHO_SUPPLIED"; cp $WHO_SUPPLIED $WHO; fi
# 2. Canada, Europe, and global results
run $S final_analysis           $CCDSS $PANEL $STATCAN $GBD_2019_2023 $CELLS $FR $ATLAS $OUT/final
run $S gbd_superregion          $LOCMAP $CELLS $OUT/final/T7_superregion_2023.csv
run $S std_certified            $CCDSS $STATCAN $OUT/display/std_certified_vs_excess.csv
run $S std_change_periods       $CCDSS $STATCAN $OUT/display/std_change_periods.json
run $S conservative_5yr         $CCDSS $STATCAN $GBD_2019_2023 $OUT/final/conservative_5yr.json
run $S conservative_ci          $CCDSS $STATCAN $GBD_2019_2023
run $S canada_yll_trend         $CCDSS $STATCAN $GBD_2019_2023 $OUT/final/canada_yll_trend.json
run $S figure_data              $CCDSS $STATCAN $GBD_2019_2023 $OUT/final $WHO $OUT/figdata
run $S table_values             $CCDSS $STATCAN $ATLAS $WHO $OUT/final/T5_international.csv $OUT/display/table_values.json
run $S display_tables           $OUT/final $CCDSS $STATCAN $GBD_2019_2023 $OUT/display
run $S last_numbers             $CCDSS $STATCAN $GBD_2019_2023 $OUT/final/T3_canada_standardised_annual.csv $OUT/display/last_numbers.json
run $S extra_numbers            $CCDSS $OUT/final/extra_numbers.json
run $S age_shift                $CCDSS $STATCAN $OUT/display/TableS4_age_share.csv
run $S canada_daly_table        $CCDSS $GBD_2019_2023 $OUT/display/canada_daly_ge40.json
run $S gbd_vs_who               $WHO $GBD_2019_2023 $OUT/final/T5_international.csv $OUT/display/gbd_vs_who_registered.csv
# 3. GBD trends, 1990-2023
run $S gbd_trends               $GBD_ANNUAL $GBD_2019_2023 $LOCMAP $WHO $OUT/display/std_certified_vs_excess.csv $OUT/final/T3_canada_standardised_annual.csv $OUT/gbd_trends
run $S gbd_annual               $GBD_ANNUAL $GBD_2019_2023 $CCDSS $OUT/gbd_trends
run $S gbd_expected             $GBD_ANNUAL $GBD_2019_2023 $LOCMAP $OUT/gbd_trends $OUT/figdata
run $S gbd_counterfactual       $GBD_ANNUAL $GBD_2019_2023 $LOCMAP $OUT/gbd_trends
run $S gbd_expected_sensitivity $GBD_ANNUAL $OUT/gbd_trends/G12_expected_sensitivity.json
run $S gbd_trend_checks         $GBD_ANNUAL $GBD_2019_2023 $STATCAN $OUT/gbd_trends
run $S gbd_trend_figdata        $OUT/gbd_trends $OUT/figdata
# 4. derived series for figures and tables (numbers only; no plotting)
run $S figure1b_population      $OUT/final/T3_canada_standardised_annual.csv $OUT/figdata/Figure1B_whole_population_index.csv
run $S derive_figure_data       $OUT $OUT/figdata
echo "done: $OUT"
