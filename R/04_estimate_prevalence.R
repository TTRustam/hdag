# -----------------------------------------------------------------------------
# Estimate empirical and smoothed SHARE prevalence
#
# Input:  data/derived/SHARE/share_prepared.rds
# Output: data/derived/SHARE/share_prevalence.rds
#
# Prevalence is estimated from observed interview states. Age is smoothed with
# a low-rank P-spline; calendar year enters linearly on the logit scale.
# Both fitting and prediction are restricted to 2007 through 2019.
# -----------------------------------------------------------------------------

if (!file.exists("R/00_functions.R")) {
  stop("Run this script from the hdag_share_msm repository root.")
}

source("R/00_functions.R")


input_file <- file.path("data", "derived", "SHARE", "share_prepared.rds")
output_dir <- file.path("data", "derived", "SHARE")
output_file <- file.path(output_dir, "share_prevalence.rds")

if (!file.exists(input_file)) {
  stop("Prepared SHARE data not found. Run R/02_prepare_SHARE.R first.")
}

share <- readRDS(input_file)
health_vars <- c("sphus", "adl", "iadl", "gali", "chronic")
year_range <- c(2007L, 2019L)

required_columns <- c("wave", "gender", "year", "age_dec", health_vars)
missing_columns <- setdiff(required_columns, names(share))

if (length(missing_columns) > 0) {
  stop(
    "Prepared SHARE data are missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

# Use the actual calendar year of interview, not the rounded decimal date.
interview_data <- share %>%
  # Death events are not prevalence observations and never enter the denominator.
  filter(wave != "death") %>%
  filter(between(year, year_range[1], year_range[2])) %>%
  filter(!is.na(gender), !is.na(year), !is.na(age_dec)) %>%
  mutate(
    sex = gender,
    age = round(age_dec)
  ) %>%
  filter(between(age, 20, 110))

if (nrow(interview_data) == 0) {
  stop("No SHARE interview observations fall in the selected years 2007 through 2019.")
}

prediction_years <- seq(year_range[1], year_range[2], by = 1)
prediction_ages <- seq(20, 110, by = 1)
sexes <- sort(unique(interview_data$sex))

empirical_parts <- vector("list", length(health_vars))
smooth_parts <- list()
model_parts <- list()
names(empirical_parts) <- health_vars
part_number <- 0L

for (health_var in health_vars) {
  outcome_data <- interview_data %>%
    transmute(
      sex,
      year,
      age,
      state = .data[[health_var]]
    ) %>%
    filter(state %in% c("H", "U"))
  
  empirical <- outcome_data %>%
    group_by(sex, year, age) %>%
    summarise(
      n_unhealthy = sum(state == "U"),
      n_healthy = sum(state == "H"),
      .groups = "drop"
    ) %>%
    mutate(
      n = n_unhealthy + n_healthy,
      prevalence = n_unhealthy / n,
      health_var = health_var,
      estimate = "empirical"
    ) %>%
    select(
      health_var,
      sex,
      year,
      age,
      estimate,
      prevalence,
      n_unhealthy,
      n_healthy,
      n
    )
  
  empirical_parts[[health_var]] <- empirical
  
  for (sex in sexes) {
    sex_data <- empirical %>%
      filter(.data$sex == .env$sex)
    
    prevalence_fit <- gam(
      prevalence ~ s(age, m = 1, k = 3, bs = "ps") + year,
      weights = n,
      family = binomial(link = "logit"),
      data = sex_data
    )
    
    prediction_data <- expand_grid(
      age = prediction_ages,
      year = prediction_years
    )
    prediction_data$prevalence <- predict(
      prevalence_fit,
      newdata = prediction_data,
      type = "response"
    )
    
    smooth <- prediction_data %>%
      mutate(
        health_var = health_var,
        sex = sex,
        estimate = "smooth",
        n_unhealthy = NA_real_,
        n_healthy = NA_real_,
        n = NA_real_
      ) %>%
      select(
        health_var,
        sex,
        year,
        age,
        estimate,
        prevalence,
        n_unhealthy,
        n_healthy,
        n
      )
    
    part_number <- part_number + 1L
    smooth_parts[[part_number]] <- smooth
    model_parts[[part_number]] <- tibble(
      health_var = health_var,
      sex = sex,
    )
  }
}

prevalence <- bind_rows(bind_rows(empirical_parts), bind_rows(smooth_parts)) %>%
  arrange(health_var, sex, year, age, estimate)
prevalence_models <- bind_rows(model_parts)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(prevalence, output_file)
saveRDS(prevalence_models, file.path(output_dir, "share_prevalence_models.rds"))

message("Empirical and smoothed prevalence saved to: ", output_file)
print(count(prevalence, health_var, sex, estimate, name = "n_rows"))
