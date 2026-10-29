# 00 load dependencies and custom functions
source("R/00_functions.R")

# 01 manage SHARE files

# SHARE data inputs:
#
# Download Release 9.0.0 SPSS files from the SHARE Data Portal:
# https://releases.sharedataportal.eu/users/login
#
# In the Releases table, download these entries (don't download wave 3):
#   "Wave 1 Release 9.0.0"  — sharew1_rel9-0-0_ALL_datasets_spss.zip
#   "Wave 2 Release 9.0.0"  — sharew2_rel9-0-0_ALL_datasets_spss.zip
#   "Wave 4 Release 9.0.0"  — sharew4_rel9-0-0_ALL_datasets_spss.zip
#   "Wave 5 Release 9.0.0"  — sharew5_rel9-0-0_ALL_datasets_spss.zip
#   "Wave 6 Release 9.0.0"  — sharew6_rel9-0-0_ALL_datasets_spss.zip
#   "Wave 7 Release 9.0.0"  — sharew7_rel9-0-0_ALL_datasets_spss.zip
#   "Wave 8 Release 9.0.0"  — sharew8_rel9-0-0_ALL_datasets_spss.zip
#   "Wave 9 Release 9.0.0"  — sharew9_rel9-0-0_ALL_datasets_spss.zip
#
# Also download this entry:
#   "All Waves Coverscreen Release 9.0.0"
#       sharewX_rel9-0-0_gv_allwaves_cv_r_spss.zip
# You should download file directly into the root of this project folder
do_this <- FALSE
if (do_this){
source("R/01_extract_SHARE_files.R")
}

# 02 recoding, reshaping, ready to estimate hazards
source("R/02_prepare_SHARE.R")

# 03 estimate hazards with msm (takes a couple minutes)
source("R/03_estimate_hazards.R")

# 04 estimate prevalence (used for adjustments)
source("R/04_estimate_prevalence.R")

# 05 read and smooth HMD mortality envelopes (used for adjustments)
# for this it's best to add your HMD credentials to your r Environment file.
# do usethis::edit_r_environ(),
# and add two entries:
# HMD_USER=your_hmd_username (notice no quotes!)
# HMD_PASSORD=your_hmd_password
# leave an empty line after, save, then restart R session.
source("R/05_smooth_HMD.R")

# 06 perform hazards adjustment using outputs of steps 3,4,5
source("R/06_adjust_hazards.R")

# 07 transform hazards to probabilities (takes a couple minutes)
source("R/07_derive_probabilities.R")

# 08 calculate h-dagger components
source("R/08_hdag.R")

# 09 calculate 2d manhattan inequality to compare with hdag
source("R/09_manhattan.R")


