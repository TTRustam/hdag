# -----------------------------------------------------------------------------
# Calculate life-table survivorship and HDAG from SHARE transition probabilities.
#
# Input:
#   data/derived/SHARE/share_probabilities_adjusted.rds
#
# Output:
#   data/derived/SHARE/share_hdag.rds
#   data/derived/SHARE/share_hdag_summary.rds
#
# The probability file contains both Unadjusted and Adjusted transitions. This
# script retains both types. mscalc calculates state survivorship (lxs) from
# the initial H/U composition and age-specific transition probabilities.
# -----------------------------------------------------------------------------

if (!file.exists("R/00_functions.R")) {
  stop("Run this script from the hdag_share_msm repository root.")
}

source("R/00_functions.R")

input_dir <- file.path("data", "derived", "SHARE")
input_file <- file.path(input_dir, "share_probabilities_adjusted.rds")
output_file <- file.path(input_dir, "share_hdag.rds")
summary_file <- file.path(input_dir, "share_hdag_summary.rds")

if (!file.exists(input_file)) {
  stop("Transition probabilities not found. Run R/07_derive_probabilities.R first: ", input_file)
}

probabilities <- readRDS(input_file)

required_columns <- c(
  "health_var", "sex", "time", "age", "from", "to", "prob", "type"
)
missing_columns <- setdiff(required_columns, names(probabilities))

if (length(missing_columns) > 0) {
  stop("Probability input is missing column(s): ", paste(missing_columns, collapse = ", "))
}

if (anyNA(probabilities$prob) ||
    any(!is.finite(probabilities$prob)) ||
    any(probabilities$prob < 0 | probabilities$prob > 1)) {
  stop("Transition probabilities must be finite values between zero and one.")
}

invalid_transitions <- probabilities %>%
  filter(
    !from %in% c("H", "U") |
      !to %in% c("H", "U", "D")
  )

if (nrow(invalid_transitions) > 0) {
  stop("Probability input must have H/U origins and H/U/D destinations.")
}

duplicate_rows <- probabilities %>%
  count(health_var, sex, time, age, from, to, type) %>%
  filter(n > 1)

if (nrow(duplicate_rows) > 0) {
  stop("Probability input has duplicate transitions within a model, sex, year, age, and type.")
}

required_transitions <- c("HH", "HU", "HD", "UH", "UU", "UD")





share_hdag <- probabilities %>%
  group_by(health_var, sex, time, type) %>%
  group_modify(
    calculate_group_hdag,
    required_transitions = required_transitions
  ) %>%
  ungroup() %>%
  arrange(health_var, sex, time, type, age)

share_hdag_summary <- share_hdag %>%
  group_by(health_var, sex, time, type) %>%
  summarise(
    across(contains("dag"), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(health_var, sex, time, type)

dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(share_hdag, output_file)
saveRDS(share_hdag_summary, summary_file)

message("Age-specific HDAG results saved to: ", output_file)
message("Summed HDAG results saved to: ", summary_file)
print(count(share_hdag, health_var, sex, type, name = "n_age_rows"))
