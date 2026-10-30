# -----------------------------------------------------------------------------
# Adjust SHARE mortality hazards to match smoothed Spanish HMD mortality rates.
#
# Inputs:
#   data/derived/SHARE/share_hazards.rds
#   data/derived/SHARE/share_prevalence.rds
#   data/derived/HMD/hmd_mx_smoothed.rds
#
# Output:
#   data/derived/SHARE/share_hazards_adjusted.rds
#
# The first-pass mortality hazard ratio is Ra = h_ud / h_hd. With smoothed
# unhealthy prevalence pi and smoothed HMD mortality rate mx, reimpute the
# mortality hazards as:
#   h_hd_adjusted = mx / ((1 - pi) + pi * Ra)
#   h_ud_adjusted = Ra * h_hd_adjusted
#
# Only H->D and U->D are replaced. The H<->U hazards remain as estimated in
# the first pass. The output retains the tidy hazard columns and adds
# type = "Adjusted" to identify these adjusted hazards.
# -----------------------------------------------------------------------------

if (!file.exists("R/00_functions.R")) {
  stop("Run this script from the hdag_share_msm repository root.")
}

source("R/00_functions.R")

input_dir <- file.path("data", "derived", "SHARE")
output_dir <- input_dir

hazard_file <- file.path(input_dir, "share_hazards.rds")
prevalence_file <- file.path(input_dir, "share_prevalence.rds")
hmd_file <- file.path("data", "derived", "HMD", "hmd_mx_smoothed.rds")
output_file <- file.path(output_dir, "share_hazards2.rds")

required_files <- c(hazard_file, prevalence_file, hmd_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Required input file(s) not found. Run the preceding pipeline steps first:\n",
    paste(missing_files, collapse = "\n")
  )
}

hazards <- readRDS(hazard_file)
prevalence <- readRDS(prevalence_file)
hmd_mx_smoothed <- readRDS(hmd_file)

hazard_keys <- c("health_var", "sex", "time", "age", "from", "to")
required_hazard_columns <- c(hazard_keys, "hazard")
required_prevalence_columns <- c(
  "health_var", "sex", "year", "age", "estimate", "prevalence"
)
required_hmd_columns <- c("sex", "time", "age", "mx_smooth")

missing_hazard_columns <- setdiff(required_hazard_columns, names(hazards))
missing_prevalence_columns <- setdiff(required_prevalence_columns, names(prevalence))
missing_hmd_columns <- setdiff(required_hmd_columns, names(hmd_mx_smoothed))

if (length(missing_hazard_columns) > 0) {
  stop("Hazard input is missing column(s): ", paste(missing_hazard_columns, collapse = ", "))
}

if (length(missing_prevalence_columns) > 0) {
  stop(
    "Prevalence input is missing column(s): ",
    paste(missing_prevalence_columns, collapse = ", ")
  )
}

if (length(missing_hmd_columns) > 0) {
  stop("HMD input is missing column(s): ", paste(missing_hmd_columns, collapse = ", "))
}

duplicate_hazards <- hazards %>%
  count(health_var, sex, time, age, from, to) %>%
  filter(n > 1)

if (nrow(duplicate_hazards) > 0) {
  stop("Hazard input has duplicate transition rows for the same model, sex, time, and age.")
}

duplicate_prevalence <- prevalence %>%
  filter(estimate == "smooth") %>%
  count(health_var, sex, year, age) %>%
  filter(n > 1)

if (nrow(duplicate_prevalence) > 0) {
  stop("Smoothed prevalence has duplicate rows for a health measure, sex, year, and age.")
}

duplicate_hmd <- hmd_mx_smoothed %>%
  count(sex, time, age) %>%
  filter(n > 1)

if (nrow(duplicate_hmd) > 0) {
  stop("Smoothed HMD mortality has duplicate rows for a sex, year, and age.")
}

# Calculate Ra from the original, first-pass SHARE mortality hazards.
first_pass_mortality <- hazards %>%
  filter(to == "D", from %in% c("H", "U")) %>%
  select(health_var, sex, time, age, from, hazard) %>%
  pivot_wider(
    names_from = from,
    values_from = hazard,
    names_prefix = "h_"
  ) |> 
  filter(between(time, 2007, 2019))

if (!all(c("h_H", "h_U") %in% names(first_pass_mortality))) {
  stop("Hazards must include first-pass H->D and U->D estimates to calculate Ra.")
}

prevalence_smooth <- prevalence %>%
  filter(estimate == "smooth") %>%
  transmute(
    health_var,
    sex,
    time = year,
    age,
    prevalence
  )

hmd_mx <- hmd_mx_smoothed %>%
  transmute(sex, time, age, mx_smooth)

adjustment_inputs <- first_pass_mortality %>%
  mutate(Ra = h_U / h_H) %>%
  select(health_var, sex, time, age, Ra) %>%
  left_join(prevalence_smooth, by = c("health_var", "sex", "time", "age")) %>%
  left_join(hmd_mx, by = c("sex", "time", "age"))

invalid_adjustment_inputs <- adjustment_inputs %>%
  filter(
    is.na(Ra) | !is.finite(Ra) | Ra < 0 |
      is.na(prevalence) | !is.finite(prevalence) |
      prevalence < 0 | prevalence > 1 |
      is.na(mx_smooth) | !is.finite(mx_smooth) | mx_smooth < 0
  )

if (nrow(invalid_adjustment_inputs) > 0) {
  stop(
    nrow(invalid_adjustment_inputs),
    " mortality-hazard rows lack valid Ra, smoothed prevalence, or smoothed HMD mx. ",
    "Check that the three inputs use the same year and age grid."
  )
}

adjustment_inputs <- adjustment_inputs %>%
  mutate(
    mortality_denominator = (1 - prevalence) + prevalence * Ra
  )

if (any(adjustment_inputs$mortality_denominator <= 0)) {
  stop("The mortality adjustment denominator must be positive for every row.")
}

adjusted_mortality <- adjustment_inputs %>%
  mutate(
    h_hd_adjusted = mx_smooth / mortality_denominator,
    h_ud_adjusted = Ra * h_hd_adjusted
  ) %>%
  select(health_var, sex, time, age, h_hd_adjusted, h_ud_adjusted) %>%
  pivot_longer(
    cols = c(h_hd_adjusted, h_ud_adjusted),
    names_to = "transition",
    values_to = "hazard"
  ) %>%
  mutate(
    from = if_else(transition == "h_hd_adjusted", "H", "U"),
    to = "D",
    type = "Adjusted"
  ) %>%
  select(health_var, sex, time, age, from, to, hazard, type)

# Keep the non-mortality transition hazards from the first-pass SHARE models.
unchanged_hazards <- hazards %>%
  filter(to != "D") %>%
  mutate(type = "Adjusted") %>%
  select(health_var, sex, time, age, from, to, hazard, type)

share_hazards_adjusted <- bind_rows(
  unchanged_hazards,
  adjusted_mortality
) %>%
  arrange(health_var, sex, time, age, from, to)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hazards |> 
  mutate(type="Unadjusted") |> 
  bind_rows(share_hazards_adjusted)|> 
  filter(between(time, 2007, 2019)) |> 
  saveRDS( output_file)


message("Adjusted and unadjusted tidy hazards saved to: ", output_file)



#   filter(time == 2018, sex == "female") |> 
#   ggplot(aes(x = age, y = hazard, linetype = type, color = interaction(from, to))) +
#   geom_line() +
#   scale_y_log10() +
#   facet_wrap(~health_var) +
#   xlim(45,110)


