# -----------------------------------------------------------------------------
# Smooth Spanish HMD mortality rates for the hazard-adjustment stage.
#
# This follows the non-bootstrap HMD smoothing model in
# incdem/R/05_mc_hmd.R (https://github.com/timriffe/incdem/blob/master/R/05_mc_hmd.R):
# Poisson deaths, an exposure offset, an age P-spline,
# and a linear calendar-year term. The fitted rates are deaths per exposure.
#
# HMD input products are downloaded with HMDHFDplus::readHMDweb() on the first
# run and cached as RDS files under data/inputs/HMD/. Set HMD_USER and
# HMD_PASSWORD in the R environment to avoid the interactive credential prompt.
#
# Output:
#   data/derived/HMD/hmd_mx_smoothed.rds
#   data/derived/HMD/hmd_mx_smoothing_models.rds
#   data/derived/HMD/hmd_mx_smoothing_check.pdf
# -----------------------------------------------------------------------------

if (!file.exists("R/00_functions.R")) {
  stop("Run this script from the hdag_share_msm repository root.")
}

source("R/00_functions.R")


input_dir <- file.path("data", "inputs", "HMD")
output_dir <- file.path("data", "derived", "HMD")

deaths_file <- file.path(input_dir, "ESP_Deaths_1x1.rds")
exposures_file <- file.path(input_dir, "ESP_Exposures_1x1.rds")
smoothed_file <- file.path(output_dir, "hmd_mx_smoothed.rds")
models_file <- file.path(output_dir, "hmd_mx_smoothing_models.rds")
diagnostic_file <- file.path(output_dir, "hmd_mx_smoothing_check.pdf")

year_range <- c(2007L, 2019L)
age_range <- c(20L, 110L)
age_interval_out <- 1L

# Keep these settings in line with the non-bootstrap incdem HMD model.
age_spline_m <- 1L
age_spline_k <- 11L

# Set TRUE to retrieve fresh products from the HMD website.
refresh_hmd_data <- FALSE

dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)



if (refresh_hmd_data || !file.exists(deaths_file)) {
  message("Reading Spanish HMD deaths (Deaths_1x1)...")
  hmd_deaths_raw <- read_hmd_product("Deaths_1x1")
  saveRDS(hmd_deaths_raw, deaths_file)
} else {
  hmd_deaths_raw <- readRDS(deaths_file)
}

if (refresh_hmd_data || !file.exists(exposures_file)) {
  message("Reading Spanish HMD exposures (Exposures_1x1)...")
  hmd_exposures_raw <- read_hmd_product("Exposures_1x1")
  saveRDS(hmd_exposures_raw, exposures_file)
} else {
  hmd_exposures_raw <- readRDS(exposures_file)
}



hmd_deaths <- prepare_hmd_product(hmd_deaths_raw, year_range, age_range) %>%
  rename(deaths = value)

hmd_exposures <- prepare_hmd_product(hmd_exposures_raw, year_range, age_range) %>%
  rename(exposure = value)

hmd_mx <- inner_join(
  hmd_deaths,
  hmd_exposures,
  by = c("year", "age", "sex")
) %>%
  filter(!is.na(deaths), !is.na(exposure), exposure > 0) %>%
  mutate(
    deaths = round(deaths, 0),
    mx_observed = deaths / exposure
  ) %>%
  arrange(sex, year, age)

if (nrow(hmd_mx) == 0) {
  stop("No usable HMD deaths and exposures found in the selected years and ages.")
}

prediction_years <- seq(year_range[1], year_range[2], by = 1L)
prediction_ages <- seq(
  age_range[1],
  age_range[2],
  by = age_interval_out
)
sexes <- sort(unique(hmd_mx$sex))

smoothed_parts <- vector("list", length(sexes))
model_parts <- vector("list", length(sexes))
names(smoothed_parts) <- sexes
names(model_parts) <- sexes

for (sex in sexes) {
  sex_data <- hmd_mx %>%
    filter(.data$sex == .env$sex)
  
  # This is the incdem model: a Poisson age spline with a linear year effect
  # and log exposure offset. Prediction exposure is one, so fitted counts are
  # fitted mortality rates (mx).
  mortality_fit <- mgcv::gam(
    deaths ~ s(age, bs = "ps", m = age_spline_m, k = age_spline_k) + year,
    offset = log(exposure),
    family = poisson(link = "log"),
    data = sex_data
  )
  
  prediction_data <- expand_grid(
    age = prediction_ages,
    year = prediction_years
  ) %>%
    mutate(exposure = 1)
  
  prediction_data$mx_smooth <- predict(
    mortality_fit,
    newdata = prediction_data,
    type = "response"
  )
  
  smoothed_parts[[sex]] <- prediction_data %>%
    transmute(
      sex = sex,
      time = year,
      age,
      mx_smooth
    )
  
  model_parts[[sex]] <- tibble(
    sex = sex,
    model = list(mortality_fit)
  )
}

hmd_mx_smoothed <- bind_rows(smoothed_parts) %>%
  left_join(
    hmd_mx %>%
      transmute(
        sex,
        time = year,
        age,
        deaths,
        exposure,
        mx_observed
      ),
    by = c("sex", "time", "age")
  ) %>%
  arrange(sex, time, age)

hmd_mx_smoothing_models <- bind_rows(model_parts)

saveRDS(hmd_mx_smoothed, smoothed_file)
saveRDS(hmd_mx_smoothing_models, models_file)


message("Smoothed HMD mortality rates saved to: ", smoothed_file)
message("HMD smoothing models saved to: ", models_file)
message("Diagnostic plot saved to: ", diagnostic_file)
print(count(hmd_mx_smoothed, sex, name = "n_estimates"))
