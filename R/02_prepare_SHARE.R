# -----------------------------------------------------------------------------
# Prepare SHARE health and interview data for the next analysis stages
#
# Inputs, created by R/01_extract_SHARE_files.R:
#   data/inputs/SHARE/health/sharew{1,2,4,5,6,7,8,9}_rel9-0-0_gv_health.sav
#   data/inputs/SHARE/sharewX_rel9-0-0_gv_allwaves_cv_r.sav
#
# This step prepares all downloaded waves except wave 3. Wave 3 is excluded
# from this analysis sequence.
#
# Output:
#   data/derived/SHARE/share_prepared.rds
# -----------------------------------------------------------------------------

# Run this script from the repository root.
if (!file.exists("R/00_functions.R")) {
  stop("Run this script from the hdag_share_msm repository root.")
}

share_dir <- file.path("data", "inputs", "SHARE")
health_dir <- file.path(share_dir, "health")
output_dir <- file.path("data", "derived", "SHARE")
output_file <- file.path(output_dir, "share_prepared.rds")

wave_ids <- c(1L, 2L, 4L:9L)
wave_labels <- paste0("w", wave_ids)

health_files <- file.path(
  health_dir,
  paste0("sharew", wave_ids, "_rel9-0-0_gv_health.sav")
)
crosswave_file <- file.path(
  share_dir,
  "sharewX_rel9-0-0_gv_allwaves_cv_r.sav"
)

missing_files <- c(health_files, crosswave_file)
missing_files <- missing_files[!file.exists(missing_files)]

if (length(missing_files) > 0) {
  stop(
    "Required SHARE input file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

source("R/00_functions.R")

# The cross-wave coverscreen supplies interview dates, age, sex, and death
# information used by gather_wave_data(). Keep Spain (country code 15).
gen_files <- read_sav(
  crosswave_file,
  col_select = c(
    mergeid,
    country,
    deceased_year,
    deceased_month,
    deceased_age,
    num_range("age_int_w", 1:9),
    num_range("int_month_w", 1:9),
    num_range("int_year_w", 1:9),
    num_range("deadoralive_w", 1:9),
    gender,
    yrbirth,
    mobirth
  )
) |>
  filter(country == 15) |>
  select(-country)

# Read the health variables from each wave's generated health module.
health_by_wave <- vector("list", length(wave_labels))
names(health_by_wave) <- wave_labels

for (i in seq_along(wave_labels)) {
  wave_label <- wave_labels[i]
  
  health_by_wave[[wave_label]] <- read_sav(
    health_files[i],
    col_select = c(
      country,
      mergeid,
      gali,
      sphus,
      adl,
      iadl,
      num_range("chronicw", 1:9)
    )
  ) |>
    filter(country == 15) |>
    select(-country)
}

# Combine each wave's health measures with interview and death information.
# gather_wave_data() recodes the health measures to H/U/D and creates dates,
# ages, and the corresponding wave-specific records.
prepared_by_wave <- vector("list", length(wave_labels))
names(prepared_by_wave) <- wave_labels

for (wave_label in wave_labels) {
  prepared_by_wave[[wave_label]] <- gather_wave_data(
    data = health_by_wave[[wave_label]],
    wave = wave_label
  )
}

# Sort each respondent's records chronologically for downstream hazard
# estimation.
share <- bind_rows(prepared_by_wave, .id = "wave") |>
  arrange(mergeid, time_int_dec)


# append deaths
n_interviews <- nrow(share)

share <- append_death_observations(
  data = share,
  gen_files = gen_files
)


dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(share, output_file)

message("Prepared SHARE data saved to: ", output_file)
print(count(share, wave, gender, name = "n_records"))

