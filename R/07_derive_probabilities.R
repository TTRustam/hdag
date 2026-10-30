# -----------------------------------------------------------------------------
# Derive adjusted transition probabilities from the adjusted SHARE hazards.
#
# Input:
#   data/derived/SHARE/share_hazards_adjusted.rds
#
# Output:
#   data/derived/SHARE/share_probabilities_adjusted.rds
#
# For each health measure, sex, year, and age, the four estimated hazards are
# placed in a 3-state generator matrix Q. Its diagonal is set to the negative
# sum of each row, and the one-year transition matrix is P = expm(Q * dt).
# Death is absorbing. Only H and U as origins are retained in the tidy output,
# as in the earlier SHARE transition-probability workflow.
# -----------------------------------------------------------------------------

if (!file.exists("R/00_functions.R")) {
  stop("Run this script from the hdag_share_msm repository root.")
}

source("R/00_functions.R")

input_dir <- file.path("data", "derived", "SHARE")
input_file <- file.path(input_dir, "share_hazards2.rds")
output_file <- file.path(input_dir, "share_probabilities_adjusted.rds")

if (!file.exists(input_file)) {
  stop("Adjusted hazards not found. Run R/06_adjust_hazards.R first.")
}

hazards <- readRDS(input_file)

required_columns <- c(
  "health_var", "sex", "time", "age", "from", "to", "hazard", "type"
)
missing_columns <- setdiff(required_columns, names(hazards))

if (length(missing_columns) > 0) {
  stop(
    "Hazard input is missing column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}


if (anyNA(hazards$hazard) ||
    any(!is.finite(hazards$hazard)) ||
    any(hazards$hazard < 0)) {
  stop("Hazards must be finite, non-missing, and non-negative.")
}

duplicate_hazards <- hazards %>%
  count(health_var, sex, time, age, from, to, type) %>%
  filter(n > 1)

if (nrow(duplicate_hazards) > 0) {
  stop("Adjusted hazard input has duplicate transition rows within a stratum.")
}

invalid_transitions <- hazards %>%
  filter(
    !from %in% c("H", "U") |
      !to %in% c("H", "U", "D") |
      from == to
  )

if (nrow(invalid_transitions) > 0) {
  stop("Hazard input contains an invalid transition. Expected H/U origins and H/U/D destinations.")
}

age_interval <- 1



share_probabilities <- hazards %>%
  group_by(health_var, sex, time, age, type) %>%
  group_modify(
    hazards_to_probabilities,
    age_interval = age_interval
  ) %>%
  ungroup() %>%
  arrange(health_var, sex, time, age, from, to)

# Check that the retained H- and U-origin probability rows sum to one.
probability_row_sums <- share_probabilities %>%
  group_by(health_var, sex, time, age, type, from) %>%
  summarise(row_sum = sum(prob), .groups = "drop")

if (any(abs(probability_row_sums$row_sum - 1) > 1e-10)) {
  stop("At least one derived H- or U-origin probability row does not sum to one.")
}

dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(share_probabilities, output_file)

message("Adjusted transition probabilities saved to: ", output_file)
print(count(share_probabilities, health_var, sex, from, to, type, name = "n_rows"))
