# -----------------------------------------------------------------------------
# Calculate the 2D healthy/unhealthy-years death distribution and Manhattan
# (rook-move) lifespan inequality from SHARE transition probabilities.
#
# Input:
#   data/derived/SHARE/share_probabilities_adjusted.rds
#
# Outputs:
#   data/derived/SHARE/share_dxh.rds
#   data/derived/SHARE/share_manhattan.rds
#
# This uses calc_dxh_fast(), with the same initial H/U composition as
# R/08_hdag.R. Add an artificial age-111 row for terminal closeout so the
# modeled transitions at age 110 are applied before everyone is counted dead.
# The 2D death distribution is aggregated over state at death.
# -----------------------------------------------------------------------------

if (!file.exists("R/00_functions.R")) {
  stop("Run this script from the hdag_share_msm repository root.")
}

source("R/00_functions.R")

input_dir <- file.path("data", "derived", "SHARE")
input_file <- file.path(input_dir, "share_probabilities_adjusted.rds")
distribution_file <- file.path(input_dir, "share_dxh.rds")
manhattan_file <- file.path(input_dir, "share_manhattan.rds")

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

radix_age <- 20L
last_age <- 110L
expected_ages <- seq(radix_age, last_age)
required_transitions <- c("HH", "HU", "HD", "UH", "UU", "UD")

probabilities_for_dxh <- probabilities %>%
  filter(between(age, radix_age, last_age))



# calc_dxh_fast() returns death mass split by state at death. Collapse that detail
# to the joint distribution d(h,u) used in the Manhattan index.
share_dxh <- probabilities_for_dxh %>%
  group_by(health_var, sex, time, type) %>%
  group_modify(
    calculate_group_dxh,
    expected_ages = expected_ages,
    required_transitions = required_transitions,
    radix_age = radix_age,
    last_age = last_age
  ) %>%
  ungroup() %>%
  group_by(health_var, sex, time, type, h, u) %>%
  summarise(
    age = first(age),
    x = first(x),
    d_hu = sum(dxsc),
    .groups = "drop"
  ) %>%
  arrange(health_var, sex, time, type, h, u)

distribution_checks <- share_dxh %>%
  group_by(health_var, sex, time, type) %>%
  summarise(total_death_mass = sum(d_hu), .groups = "drop")

if (any(abs(distribution_checks$total_death_mass - 1) > 1e-7)) {
  stop("At least one 2D death distribution does not sum to one after terminal-age closeout.")
}

distribution_means <- share_dxh %>%
  group_by(health_var, sex, time, type) %>%
  summarise(
    mean_h = sum(h * d_hu),
    mean_u = sum(u * d_hu),
    mean_lifespan = sum(x * d_hu),
    .groups = "drop"
  )

share_manhattan <- share_dxh %>%
  left_join(distribution_means, by = c("health_var", "sex", "time", "type")) %>%
  group_by(health_var, sex, time, type) %>%
  summarise(
    mean_h = first(mean_h),
    mean_u = first(mean_u),
    mean_lifespan = first(mean_lifespan),
    manhattan_h = sum(abs(h - mean_h) * d_hu),
    manhattan_u = sum(abs(u - mean_u) * d_hu),
    .groups = "drop"
  ) %>%
  mutate(manhattan = manhattan_h + manhattan_u) %>%
  arrange(health_var, sex, time, type)

dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(share_dxh, distribution_file)
saveRDS(share_manhattan, manhattan_file)

message("2D healthy/unhealthy-years death distributions saved to: ", distribution_file)
message("Manhattan inequality and H/U components saved to: ", manhattan_file)
print(share_manhattan)

