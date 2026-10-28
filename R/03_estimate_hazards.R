# -----------------------------------------------------------------------------
# Estimate SHARE multi-state transition hazards and save them in tidy form
#
# Input:  data/derived/SHARE/share_prepared.rds
# Output: data/derived/SHARE/share_hazards.rds
#         data/derived/SHARE/share_msm_models.rds
#
# Calendar time is modeled as a linear covariate for this first fit.
# msm holds a time-varying covariate fixed between observations, so a simple
# 2020/2021 indicator would not represent within-gap COVID exposure in extended
# wave-9 intervals. Revisit that specification after checking interview dates.
# -----------------------------------------------------------------------------

if (!file.exists("R/00_functions.R")) {
  stop("Run this script from the hdag_share_msm repository root.")
}

source("R/00_functions.R")

input_file <- file.path("data", "derived", "SHARE", "share_prepared.rds")
output_dir <- file.path("data", "derived", "SHARE")
hazard_file <- file.path(output_dir, "share_hazards.rds")
model_file <- file.path(output_dir, "share_msm_models.rds")

if (!file.exists(input_file)) {
  stop("Prepared SHARE data not found. Run R/02_prepare_SHARE.R first.")
}

share <- readRDS(input_file)

health_vars <- c("sphus", "adl", "iadl", "gali", "chronic")

required_columns <- c(
  "wave", "mergeid", "gender", "age_dec", "time_int_dec", health_vars
)
missing_columns <- setdiff(required_columns, names(share))

if (length(missing_columns) > 0) {
  stop(
    "Prepared SHARE data are missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

if (!any(share$wave == "death", na.rm = TRUE)) {
  stop(
    "No appended death rows found. Rebuild the prepared data by running ",
    "R/02_prepare_SHARE.R first."
  )
}

# Healthy and unhealthy are transient states; death is absorbing.
Q <- rbind(
  c(0,   0.1, 0.1),
  c(0.1, 0,   0.1),
  c(0,   0,   0)
)
diag(Q) <- -rowSums(Q)

model_parts <- vector("list", length(health_vars))
names(model_parts) <- health_vars

for (health_var in health_vars) {
  model_parts[[health_var]] <- fit_msm_health(
    data = share,
    health_var = health_var,
    Q = Q
  )
}

models <- bind_rows(model_parts)

if (nrow(models) == 0) {
  stop("No SHARE MSMs were fitted. Check the prepared data and health states.")
}

# Predict instantaneous transition hazards at integer ages and calendar years
# spanning the observed follow-up. These are generator intensities, not
# transition probabilities; the later adjustment step will modify hazards.
prediction_ages <- seq(20, 110, by = 1)
prediction_years <- seq(
  floor(min(share$time_int_dec, na.rm = TRUE)),
  ceiling(max(share$time_int_dec, na.rm = TRUE)),
  by = 1
)

# Each fitted model contributes four allowed transition hazards: H->U, H->D,
# U->H, and U->D. Death is absorbing, so there are no outgoing D transitions.
hazard_parts <- list()
part_number <- 0L

for (model_row in seq_len(nrow(models))) {
  model_fit <- models$model[[model_row]]
  health_var <- models$health_var[[model_row]]
  sex <- models$gender[[model_row]]
  
  for (prediction_year in prediction_years) {
    for (prediction_age in prediction_ages) {
      q_matrix <- qmatrix.msm(
        x = model_fit,
        covariates = list(
          age_dec = prediction_age,
          time_int_dec = prediction_year
        ),
        ci = "none"
      )
      
      q_tidy <- as.data.frame(as.table(q_matrix)) %>%
        rename(from = Var1, to = Var2, hazard = Freq) %>%
        mutate(
          from = recode(
            as.character(from),
            "State 1" = "H",
            "State 2" = "U",
            "State 3" = "D"
          ),
          to = recode(
            as.character(to),
            "State 1" = "H",
            "State 2" = "U",
            "State 3" = "D"
          )
        ) %>%
        filter(from != to, hazard > 0) %>%
        transmute(
          health_var = health_var,
          sex = sex,
          time = prediction_year,
          age = prediction_age,
          from,
          to,
          hazard = as.numeric(hazard)
        )
      
      part_number <- part_number + 1L
      hazard_parts[[part_number]] <- q_tidy
    }
  }
}

hazards <- bind_rows(hazard_parts) %>%
  arrange(health_var, sex, time, age, from, to)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(models, model_file)
saveRDS(hazards, hazard_file)

message("Fitted MSMs saved to: ", model_file)
message("Tidy transition hazards saved to: ", hazard_file)
print(count(hazards, health_var, sex, from, to, name = "n_estimates"))
