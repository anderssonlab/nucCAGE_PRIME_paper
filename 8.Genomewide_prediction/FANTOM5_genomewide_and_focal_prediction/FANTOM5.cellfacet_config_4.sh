# config.sh

### 1 ### _1_get_ctss_from_bw.r
CAGE_BW_DIR="/projects/ralab/data/projects/nucleiCAGEproject/0.External_resources/FANTOM5/pooled_bw"
DESIGN_MATRIX="FANTOM5.cellfacet_design_matrix_4.tsv"
OUTPUT_DIR="FANTOM5.cellfacet_4"
CTSS_RSE_NAME="ctss_rse_FANTOM5.cellfacet.pooled.bw_4.rds"

### 0 ### _0_validate_ctss_and_region.r
CTSS_RSE_RDS="/FANTOM5.cellfacet_4/ctss_rse_FANTOM5.cellfacet.pooled.bw_4.rds"
REGION_RDS="/FANTOM5.rmSingletons/FANTOM5_PRIME_pred_all_combined_thresh0.75_d0.1.rds"

### 4 ### _4_get_profile.r
# CTSS_RSE_RDS, REGION_RDS, OUTPUT_DIR from above
PROFILE_MAIN_DIR="FANTOM5_focal_profiles"
PROFILE_FORMAT="npz"
# add --save_count_profiles if you want to save count profiles

### 5 ### _5_predict_profile_probability.py
# OUTPUT_DIR, PROFILE_MAIN_DIR, and PROFILE_FORMAT from above
PYTHON_PATH="/usr/bin/python3"
MODEL_PATH=$(Rscript -e 'cat(system.file("model", "PRIME_GM12878_model_1.0.sav", package = "PRIMEmodel"))')
PREFIX_OUT_NAME="FANTOM5-cellfacet-on-PRIMEmodel"
