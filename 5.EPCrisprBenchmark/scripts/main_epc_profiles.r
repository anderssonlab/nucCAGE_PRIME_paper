# Working directory and environment variables
work_dir <- this.path::this.dir()
setwd(work_dir)
dotenv::load_dot_env("epcrispr_benchmark.env")

# Source additional script
source("../../CAGEfightR_extensions/heatmap.R")
source("../r_utils/utils.r")
source("../r_gr_management/manipulate_gr.r")
source("../r_profiles/directories.r")
source("../r_profiles/profiles.r")



# 0. setting

# Retrieve environment variables
ext_dis <- as.numeric(Sys.getenv("EXT_DIS"))
len_vec <- ext_dis * 2 + 1

# main directories
dir_resources <- Sys.getenv("DIR_RESOURCES")
dir_results <- Sys.getenv("DIR_RESULTS")

# epc profiles directory's names
outdir_dir <- Sys.getenv("OUTPUT_DIR_EPCPROFILES")
outdir_dir_name <- Sys.getenv("OUTPUT_DIR_NAME_EPCPROFILES")

# set this according to the data
outdir_main_name <- c("metadata",
                      "metadata_subtnorm",
                      "profiles",
                      "profiles_subtnorm")
outdir_subdir_name <- c("pos", "neg")
outfile_pos_gr <- Sys.getenv("OUTFILE_POS_GR")
outfile_neg_notsig_gr <- Sys.getenv("OUTFILE_NEG_NOTSIG_GR")

# read in ranges data
pos_gr <- readRDS(file.path(dir_results, outfile_pos_gr))
neg_gr <- readRDS(file.path(dir_results, outfile_neg_notsig_gr))

# read in CTSSs data
infile_ctss_rse <- Sys.getenv("INFILE_CTSS_RSE")
ctss_rse <- readRDS(file.path(dir_resources, infile_ctss_rse))

# 1. create the output dir and subfolder
prep_profile_dir(output_dir = outdir_dir,
                 output_dir_name = outdir_dir_name,
                 output_main_name = outdir_main_name,
                 output_subdir_name = outdir_subdir_name)

# 2. create profile
# modify region_gr and outdir_subdir_name
report_time_execution(wrapup_make_profiles(ctss_rse, pos_gr,
                                           dir_results,
                                           outdir_dir_name,
                                           outdir_subdir_name[1],
                                           ext_dis, len_vec))
report_time_execution(wrapup_make_profiles(ctss_rse, neg_gr,
                                           dir_results,
                                           outdir_dir_name,
                                           outdir_subdir_name[2],
                                           ext_dis, len_vec))
