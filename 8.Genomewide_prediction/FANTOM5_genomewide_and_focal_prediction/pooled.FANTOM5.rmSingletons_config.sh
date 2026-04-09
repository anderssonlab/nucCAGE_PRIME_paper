# config.sh

### 1 ### _1_get_ctss_from_bw.r
CAGE_BW_DIR="/maps/projects/ralab/data/projects/nucleiCAGEproject/0.External_resources/FANTOM5/"
DESIGN_MATRIX="pooled.FANTOM5.rmSingletons_design_matrix.tsv"
OUTPUT_DIR="FANTOM5.rmSingletons"
CTSS_RSE_NAME="ctss_rse_FANTOM5.rmSingletons.rds"
# add -k for keeping only standard chromosomes

### 2 ### _2_get_tc_from_ctss.r
# EXTENSION_DISTANCE=200 was fixed in the script
# OUTPUT_DIR and CTSS_RSE_NAME from above
TC_GRL_NAME="tc_grl_FANTOM5.rmSingletons.rds"

### 3 ### _3_get_sld_window_from_tc.r
# OUTPUT_DIR and TC_GRL_NAME from above
SLD_TC_GRL_NAME="tc_grl_sld_FANTOM5.rmSingletons.rds"
SLD_WINDOW=20

### 4 ### _4_get_profile.r
# CTSS_RSE_NAME, TC_GRL_NAME, OUTPUT_DIR from above
PROFILE_MAIN_DIR="FANTOM5.rmSingletons"
PROFILE_FILE_TYPE="npz"
# add --save_count_profiles if you want to save count profiles

### 5 ### _5_predict_profile_probability.py
# OUTPUT_DIR, PROFILE_MAIN_DIR, PROFILE_SUB_DIR, and PROFILE_FILE_TYPE from above
PYTHON_PATH="/usr/bin/python3"
MODEL_PATH=$(Rscript -e 'cat(system.file("model", "PRIME_GM12878_model_1.0.sav", package = "PRIMEmodel"))')
PREFIX_OUT_NAME="FANTOM5_PRIME"


## 6 ### _6_apply_post_processing_coreovlwith-d.r
# OUTPUT_DIR from above
PARTIAL_NAME="pred_all"
THRESHOLD=0.75
SCORE_DIFF=0.10
CORE_WIDTH=150
