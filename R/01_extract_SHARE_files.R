# -----------------------------------------------------------------------------
# SHARE data inputs
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
#
# This workflow uses the gv_health module from each wave bundle and the
# gv_allwaves_cv_r module from the All Waves Coverscreen bundle. Wave 3 is
# excluded from the current analysis sequence.
#
# After extraction, the required .sav files are placed in:
#   data/inputs/SHARE/health/   (wave-specific gv_health files)
#   data/inputs/SHARE/          (cross-wave gv_allwaves_cv_r file)
#
# SHARE data are not committed to this repository. Keep downloaded ZIPs and
# extracted data out of version control.
# -----------------------------------------------------------------------------

if (!file.exists("R/00_functions.R")) {
  stop("Run this script from the hdag_share_msm repository root.")
}

source("R/00_functions.R")

download_dir <- getwd() # change to path.expand("~/Downloads") if you downloaded there...
project_dir <- getwd()

if (!dir.exists(download_dir)) {
  stop("Download folder not found: ", download_dir)
}

share_dir <- file.path(project_dir, "data", "inputs", "SHARE")
health_dir <- file.path(share_dir, "health")

dir.create(health_dir, recursive = TRUE, showWarnings = FALSE)

health_files <- sprintf(
  "sharew%d_rel9-0-0_gv_health.sav",
  c(1, 2, 4:9)
)

crosswave_file <- "sharewX_rel9-0-0_gv_allwaves_cv_r.sav"
targets <- c(health_files, crosswave_file)

zip_files <- list.files(
  download_dir,
  pattern = "\\.zip$",
  recursive = TRUE,
  full.names = TRUE
)

sources <- list()

for (target in targets) {
  source <- find_share_file(target, download_dir, zip_files)
  
  if (is.null(source)) {
    stop("Could not find ", target, " in the download folder or its ZIP files.")
  }
  
  sources[[target]] <- source
}

for (target in targets) {
  source <- sources[[target]]
  
  if (target %in% health_files) {
    destination <- file.path(health_dir, target)
  } else {
    destination <- file.path(share_dir, target)
  }
  
  if (source$type == "file") {
    copied <- file.copy(source$path, destination, overwrite = TRUE)
    
    if (!copied) {
      stop("Could not copy ", target, " to ", destination)
    }
  } else {
    staging_dir <- tempfile("share_extract_")
    dir.create(staging_dir)
    
    utils::unzip(
      source$path,
      files = source$member,
      exdir = staging_dir,
      junkpaths = TRUE,
      overwrite = TRUE
    )
    
    staged_file <- file.path(staging_dir, basename(source$member))
    
    if (!file.exists(staged_file)) {
      unlink(staging_dir, recursive = TRUE)
      stop("Could not extract ", target, " from ", basename(source$path))
    }
    
    copied <- file.copy(staged_file, destination, overwrite = TRUE)
    unlink(staging_dir, recursive = TRUE)
    
    if (!copied) {
      stop("Could not copy ", target, " to ", destination)
    }
  }
  
  message("Ready: ", destination)
}

health_ready <- list.files(
  health_dir,
  pattern = "\\.sav$",
  full.names = FALSE
)

if (!setequal(health_ready, health_files)) {
  stop("The health folder does not contain exactly the eight expected files.")
}

if (!file.exists(file.path(share_dir, crosswave_file))) {
  stop("The cross-wave file is missing from ", share_dir)
}

message("SHARE inputs are ready in: ", share_dir)
