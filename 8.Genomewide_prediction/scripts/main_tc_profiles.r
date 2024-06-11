# Working directory and environment variables
work_dir <- this.path::this.dir()
setwd(work_dir)
dotenv::load_dot_env("genomewide.env")

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

# profiles directory's names
outdir_dir <- Sys.getenv("OUTPUT_DIR_PROFILES")
outdir_dir_name <- Sys.getenv("OUTPUT_DIR_NAME_PROFILES")

# set this according to the data
outdir_main_name <- c("metadata",
                      "metadata_subtnorm",
                      "profiles",
                      "profiles_subtnorm")
outdir_subdir_name <- c("tcs")

# read in ranges data
outfile_tc_grl <- Sys.getenv("OUTFILE_TC_GRL")
tc_grl <- readRDS(file.path(dir_results, outfile_tc_grl))

# read in CTSSs data
infile_ctss_rse <- Sys.getenv("INFILE_CTSS_RSE")
ctss_rse <- readRDS(file.path(dir_resources, infile_ctss_rse))
# if you wanna choose some libraries, you can do it here:
ctss_rse <- ctss_rse[, c(7, 8, 13, 14, 17)] # nolint



# 1. create the output dir and subfolder

prep_profile_dir(output_dir = outdir_dir,
                 output_dir_name = outdir_dir_name,
                 output_main_name = outdir_main_name,
                 output_subdir_name = outdir_subdir_name)



# 2. create profile

# modify region_gr and outdir_subdir_name
report_time_execution(wrapup_make_profiles(ctss_rse, tc_grl,
                                           dir_results,
                                           outdir_dir_name,
                                           outdir_subdir_name[1],
                                           ext_dis, len_vec))
