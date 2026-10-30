# -----------------------------------------------------------------------------
# Shared functions and package setup for the SHARE H/DAG workflow
#
# Manifest: function, purpose, and calling script
#   gather_wave_data()           Recode one wave's SHARE interview health data; 02_prepare_SHARE.R
#   append_death_observations()  Append observed terminal deaths; 02_prepare_SHARE.R
#   fit_msm_health()             Fit sex-specific MSMs for one health measure; 03_estimate_hazards.R
#   init_constant()              Set initial H/U composition; 08_hdag.R, 09_manhattan.R
#   calc_dxh_fast()              Build the 2D healthy/unhealthy death mass; 09_manhattan.R
#   find_share_file()            Find a SHARE module in files or ZIPs; 01_extract_SHARE_files.R
#   read_hmd_product()           Read a Spanish HMD product; 05_smooth_HMD.R
#   prepare_hmd_product()        Tidy HMD deaths or exposures; 05_smooth_HMD.R
#   hazards_to_probabilities()   Convert one hazard stratum to probabilities; 07_derive_probabilities.R
#   calculate_group_hdag()       Calculate life table and HDAG; 08_hdag.R
#   calculate_group_dxh()        Calculate one stratum's 2D death distribution; 09_manhattan.R
#
#   Legacy helpers not used by steps 01-09 remain commented out at the end.
#
# Package setup: install missing CRAN packages and timriffe/mscalc from GitHub,
# then attach the packages used across the numbered scripts. This requires an
# internet connection on a fresh machine; pak is installed from CRAN if needed.
# -----------------------------------------------------------------------------

cran_packages <- c(
  "tidyverse",
  "haven",
  "msm",
  "mgcv",
  "expm",
  "janitor",
  "HMDHFDplus"
)

github_packages <- c(
  mscalc = "github::timriffe/mscalc"
)

if (!requireNamespace("pak", quietly = TRUE)) {
  utils::install.packages("pak", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("pak", quietly = TRUE)) {
  stop("Could not install 'pak'. Install it from CRAN, then source this file again.")
}

missing_cran_packages <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
missing_github_packages <- github_packages[
  !vapply(names(github_packages), requireNamespace, logical(1), quietly = TRUE)
]
packages_to_install <- c(missing_cran_packages, unname(missing_github_packages))

if (length(packages_to_install) > 0) {
  pak::pkg_install(packages_to_install)
}

packages_to_load <- c(cran_packages, names(github_packages))
for (package_name in packages_to_load) {
  suppressPackageStartupMessages(
    library(package_name, character.only = TRUE)
  )
}



#' Gather and recode one SHARE wave's interview health records.
#' Used by: R/02_prepare_SHARE.R.
gather_wave_data <- function(data = health_data[[6]],
                             wave = "w7") {
  
  disease <- str_c("chronic",       wave)
  age_int <- str_c("age_int_",      wave)
  int_m   <- str_c("int_month_",    wave)
  int_y   <- str_c("int_year_",     wave)
  wave    <- str_c("deadoralive_" , wave)
  
  sample <- gen_files %>%
    dplyr::select(mergeid,
                  gender,
                  all_of(c(int_m, int_y, age_int)),
                  contains(wave),
                  yrbirth,
                  mobirth) %>%
    filter(!!sym(int_y) > 0) %>%
    set_names("mergeid", 
              "gender", 
              "int_month",
              "int_year",
              "age",
              "deadoralive",
              "yrbirth",
              "mobirth") %>%
    mutate(int_month = ifelse(int_month < 0 | is.na(int_month), 1, int_month)) %>%
    filter(deadoralive %in% c(1, 2)) %>% 
    dplyr::select(-deadoralive) %>%
    unite(time_int, c("int_year", "int_month"), sep = ".") %>% 
    mutate(time_int = ym(time_int))
  
  generated <- sample %>%
    filter(age > 0) %>%
    mutate(
      age         = as.numeric(age),
      gender      = ifelse(gender == 1, "male", "female"),
      year        = year(time_int),
    ) %>%
    unite(time_born, c("yrbirth", "mobirth"), sep = ".") %>%
    mutate(time_born = ym(time_born)) %>%
    mutate(time_born_dec = decimal_date(time_born),
           time_int_dec  = decimal_date(time_int),
           age_dec       = time_int_dec - time_born_dec) %>%
    select(-c(time_int, time_born))
  
  # age_years = interval(time_born, time_int) / years(1)
  
  health <- data %>%
    rename("chronic" := !!sym(disease))
  
  health %>%
    left_join(generated, by = join_by(mergeid)) %>% 
    # prepare 3 health definitions
    mutate(gali = case_when(
      gali == 0 ~ "H",
      gali == 1 ~ "U",
      TRUE      ~ NA_character_
    )) %>%
    mutate(sphus = case_when(
      sphus %in% c(1, 2, 3) ~ "H",
      sphus %in% c(4, 5)    ~ "U",
      TRUE                  ~ NA_character_
    )) %>%
    mutate(chronic = case_when(
      chronic <  2  ~ "H",
      chronic >= 2  ~ "U",
      TRUE  ~ NA_character_
    )) %>%
    mutate(adl = case_when(
      adl == 0 ~ "H",
      adl  > 0 ~ "U",
      TRUE ~ NA_character_
    )) %>%
    mutate(iadl = case_when(
      iadl == 0 ~ "H",
      iadl  > 0 ~ "U",
      TRUE ~ NA_character_
    )) %>%
    # remove cases where none of the health states is known
    # Maybe fill??
    filter(!if_all(c(chronic, gali, sphus, adl, iadl), is.na))
  
}

#' Append valid death records after each person's last interview.
#' Used by: R/02_prepare_SHARE.R.
append_death_observations <- function(data, gen_files) {
  
  last_interview <- data %>%
    filter(!is.na(time_int_dec)) %>%
    group_by(mergeid) %>%
    summarise(last_interview_dec = max(time_int_dec), .groups = "drop")
  
  death_observations <- gen_files %>%
    transmute(
      mergeid,
      gender_code = as.numeric(gender),
      birth_year = as.numeric(yrbirth),
      birth_month = as.numeric(mobirth),
      death_year = as.numeric(deceased_year),
      death_month = as.numeric(deceased_month),
      death_age = as.numeric(deceased_age),
      death_status = as.numeric(deadoralive_w9)
    ) %>%
    filter(
      death_status == 2,
      death_year > 0,
      death_age >= 0
    ) %>%
    mutate(
      death_month = ifelse(
        is.na(death_month) | death_month < 1 | death_month > 12,
        1,
        death_month
      ),
      gender = ifelse(gender_code == 1, "male", "female"),
      time_dead = ym(sprintf("%04d-%02d", death_year, death_month)),
      time_born = ym(if_else(
        !is.na(birth_year) & birth_year > 0 & birth_month %in% 1:12,
        sprintf("%04d-%02d", birth_year, birth_month),
        NA_character_
      )),
      time_int_dec = decimal_date(time_dead),
      time_born_dec = decimal_date(time_born),
      age_dec = time_int_dec - time_born_dec,
      year = year(time_dead)
    ) %>%
    inner_join(last_interview, by = "mergeid") %>%
    filter(time_int_dec > last_interview_dec) %>%
    transmute(
      wave = "death",
      mergeid,
      gender,
      age = death_age,
      year,
      time_born_dec,
      time_int_dec,
      age_dec,
      chronic = "D",
      gali = "D",
      sphus = "D",
      adl = "D",
      iadl = "D"
    )
  
  bind_rows(data, death_observations) %>%
    arrange(mergeid, time_int_dec)
}

#' Fit sex-specific continuous-time MSMs for one health measure.
#' Used by: R/03_estimate_hazards.R.
fit_msm_health <- function(data = share, 
                           health_var = "sphus", Q,
                           process_time = "age_dec",
                           pci = NULL) {
  
  model_data <- data %>%
    transmute(
      mergeid,
      gender,
      age_dec,
      time_int_dec,
      health = .data[[health_var]]
    ) %>%
    filter(
      !is.na(health),
      !is.na(time_int_dec),
      health %in% c("H", "U", "D")
    ) %>%
    arrange(mergeid, time_int_dec) %>%
    group_by(mergeid) %>%
    mutate(first_death = match("D", health)) %>%
    # Keep the first death row and discard anything after death.
    filter(is.na(first_death) | row_number() <= first_death) %>%
    ungroup() %>%
    select(-first_death) %>%
    # msm needs at least two usable observations per respondent.
    filter(n() > 1, .by = mergeid) %>%
    arrange(mergeid, time_int_dec) %>%
    ungroup()
  
  missing_age_ids <- model_data %>%
    filter(any(is.na(age_dec)), .by = mergeid) %>%
    distinct(mergeid)
  
  if (nrow(missing_age_ids) > 0) {
    warning(
      health_var,
      ": excluding ",
      nrow(missing_age_ids),
      " respondent(s) with multiple observations but missing age at one or more observations.",
      call. = FALSE
    )
  }
  
  model_data <- model_data %>%
    # Keep complete respondent histories; otherwise a missing-age death row
    # could silently turn a known death into apparent censoring.
    filter(all(!is.na(age_dec)), .by = mergeid) %>%
    arrange(mergeid, time_int_dec) %>%
    mutate(
      health = case_when(
        health == "H" ~ 1L,
        health == "U" ~ 2L,
        health == "D" ~ 3L
      ),
      # Interviews are panel snapshots; appended deaths use the reported
      # death month (the preparation helper imputes January when month is missing).
      obstype = if_else(health == 3L, 3L, 1L)
    )
  
  sexes <- sort(unique(model_data$gender))
  model_rows <- vector("list", length(sexes))
  names(model_rows) <- sexes
  model_formula <- reformulate(process_time, response = "health")
  
  for (sex in sexes) {
    sex_data <- model_data %>%
      filter(gender == sex)
    
    model_fit <- msm(
      formula    = model_formula,
      subject    = mergeid,
      data       = sex_data,
      qmatrix    = Q,
      obstype    = obstype,
      covariates = ~ age_dec + time_int_dec,
      pci        = pci,
      control    = list(
        fnscale = 5000,
        maxit = 25000
      ),
      method = "BFGS",
      gen.inits = TRUE
    )
    
    model_rows[[sex]] <- tibble(
      gender = sex,
      data = list(sex_data),
      model = list(model_fit),
      health_var = health_var
    )
  }
  
  bind_rows(model_rows)
}

#' Calculate the starting H/U mix from first-age transition probabilities.
#' Used by: calculate_group_hdag() in R/08_hdag.R and
#' calculate_group_dxh() in R/09_manhattan.R.
init_constant <- function(x) {
  transition_values <- unlist(x[c("HH", "UH", "HU", "UU")])
  transient_matrix <- matrix(transition_values, nrow = 2, byrow = TRUE)
  eigen_result <- eigen(transient_matrix)
  dominant_index <- which.max(Re(eigen_result$values))
  initial <- Re(eigen_result$vectors[, dominant_index])

  if (sum(initial) < 0) {
    initial <- -initial
  }
  if (any(initial < -1e-8) || sum(initial) <= 0) {
    stop("Could not calculate a valid initial H/U composition from the first-age probabilities.")
  }

  initial[initial < 0] <- 0
  initial <- initial / sum(initial)
  setNames(as.numeric(initial), c("H", "U"))
}

#' Calculate the 2D distribution of healthy and unhealthy years at death.
#' Used by: calculate_group_dxh() in R/09_manhattan.R.
calc_dxh_fast <- function(p_tibble,
                          age_int = 1,
                          init = c(H = 1, U = 0),
                          keep_zero = FALSE) {
  p_tibble <- p_tibble |>
    dplyr::arrange(age)
  
  ages <- p_tibble$age
  min_age <- min(ages)
  n_ages <- length(ages)
  
  H <- unname(init["H"])
  U <- unname(init["U"])
  out <- vector("list", n_ages)
  
  for (i in seq_len(n_ages)) {
    age_i <- ages[i]
    x_i <- age_i - min_age
    h_i <- seq.int(0, length(H) - 1L) * age_int
    u_i <- x_i - h_i
    terminal <- i == n_ages
    
    qH <- if (terminal) 1 else p_tibble$HD[i]
    qU <- if (terminal) 1 else p_tibble$UD[i]
    
    out[[i]] <- dplyr::bind_rows(
      tibble::tibble(
        current_state = "H",
        age = age_i,
        x = x_i,
        h = h_i,
        u = u_i,
        lxsc = H,
        dxsc = H * qH
      ),
      tibble::tibble(
        current_state = "U",
        age = age_i,
        x = x_i,
        h = h_i,
        u = u_i,
        lxsc = U,
        dxsc = U * qU
      )
    )
    
    if (!terminal) {
      # Arrival in H increments accumulated healthy duration by one step.
      H_next <- numeric(length(H) + 1L)
      H_next[-1L] <- H * p_tibble$HH[i] + U * p_tibble$UH[i]
      
      # Arrival in U leaves accumulated healthy duration unchanged.
      U_next <- numeric(length(U) + 1L)
      U_next[seq_along(U)] <- U * p_tibble$UU[i] + H * p_tibble$HU[i]
      
      H <- H_next
      U <- U_next
    }
  }
  
  ans <- dplyr::bind_rows(out)
  if (!keep_zero) {
    ans <- dplyr::filter(ans, lxsc > 0 | dxsc > 0)
  }
  ans
}

#' Find a requested SHARE data file loose or inside a ZIP archive.
#' Used by: R/01_extract_SHARE_files.R.
find_share_file <- function(target, download_dir, zip_files) {
  loose_files <- list.files(
    download_dir,
    recursive = TRUE,
    full.names = TRUE
  )
  loose_match <- loose_files[basename(loose_files) == target]
  
  if (length(loose_match) > 0) {
    return(list(type = "file", path = loose_match[1]))
  }
  
  for (zip_path in zip_files) {
    archive_names <- utils::unzip(zip_path, list = TRUE)$Name
    archive_match <- archive_names[basename(archive_names) == target]
    
    if (length(archive_match) > 0) {
      return(list(
        type = "zip",
        path = zip_path,
        member = archive_match[1]
      ))
    }
  }
  
  NULL
}

#' Read one Spanish HMD data product, using environment credentials when set.
#' Used by: R/05_smooth_HMD.R.
read_hmd_product <- function(item,
                             hmd_user = Sys.getenv("HMD_USER"),
                             hmd_password = Sys.getenv("HMD_PASSWORD")) {
  if (nzchar(hmd_user) && nzchar(hmd_password)) {
    HMDHFDplus::readHMDweb(
      CNTRY = "ESP",
      item = item,
      username = hmd_user,
      password = hmd_password
    )
  } else {
    HMDHFDplus::readHMDweb(CNTRY = "ESP", item = item)
  }
}

#' Convert one HMD product to tidy year-age-sex-value records.
#' Used by: R/05_smooth_HMD.R.
prepare_hmd_product <- function(data, year_range, age_range) {
  data %>%
    janitor::clean_names() %>%
    transmute(
      year = as.integer(readr::parse_number(as.character(year))),
      age = readr::parse_number(as.character(age)),
      male = readr::parse_number(as.character(male)),
      female = readr::parse_number(as.character(female))
    ) %>%
    filter(
      between(year, year_range[1], year_range[2]),
      between(age, age_range[1], age_range[2])
    ) %>%
    pivot_longer(
      cols = c(male, female),
      names_to = "sex",
      values_to = "value"
    )
}

#' Convert one grouped set of transition hazards to probabilities.
#' Used by: R/07_derive_probabilities.R.
hazards_to_probabilities <- function(.x, .y, age_interval) {
  states <- c("H", "U", "D")
  
  Q <- matrix(
    0,
    nrow = length(states),
    ncol = length(states),
    dimnames = list(from = states, to = states)
  )
  
  for (transition_row in seq_len(nrow(.x))) {
    from_state <- as.character(.x$from[transition_row])
    to_state <- as.character(.x$to[transition_row])
    Q[from_state, to_state] <- .x$hazard[transition_row]
  }
  
  # Death has no outgoing hazards, so its row remains zero and P_DD = 1.
  diag(Q) <- -rowSums(Q)
  
  P <- expm::expm(Q * age_interval)
  
  as.data.frame(as.table(P)) %>%
    rename(from = Var1, to = Var2, prob = Freq) %>%
    filter(from != "D") %>%
    mutate(prob = as.numeric(prob))
}

#' Calculate survivorship, state expectancy, and HDAG for one stratum.
#' Used by: R/08_hdag.R.
calculate_group_hdag <- function(.x, .y, required_transitions) {
  transition_probabilities <- .x %>%
    mutate(from_to = paste0(as.character(from), as.character(to))) %>%
    select(age, from_to, prob)
  
  probability_wide <- transition_probabilities %>%
    pivot_wider(names_from = from_to, values_from = prob) %>%
    arrange(age)
  
  missing_transitions <- setdiff(required_transitions, names(probability_wide))
  if (length(missing_transitions) > 0) {
    stop(
      "A model group is missing transition(s): ",
      paste(missing_transitions, collapse = ", ")
    )
  }
  
  if (anyNA(probability_wide[required_transitions])) {
    stop("A model group has an incomplete age-by-transition probability grid.")
  }
  
  age_steps <- diff(probability_wide$age)
  if (length(age_steps) > 0 && any(abs(age_steps - 1) > 1e-8)) {
    stop("A model group must have consecutive one-year ages for this life-table calculation.")
  }
  
  initial <- init_constant(probability_wide[1, ])
  
  # mscalc returns start-of-interval survivorship for each transient state.
  # The probabilities at the final age are retained for final-age transfers;
  # they do not change survivorship at that same age.
  lxs <- mscalc:::calc_survivorship(
    transitions = transition_probabilities,
    init = initial,
    delim = "",
    age_col = "age",
    trans_col = "from_to",
    p_col = "prob"
  ) %>%
    select(age, H, U) %>%
    rename(lh = H, lu = U)
  
  # Expected future years in H and U conditional on each starting state.
  # This backward recursion is the age-structured equivalent of the
  # fundamental-matrix calculation in the original HDAG script.
  n_ages <- nrow(probability_wide)
  expected_years <- matrix(
    0,
    nrow = n_ages,
    ncol = 4,
    dimnames = list(NULL, c("HLEh", "ULEh", "HLEu", "ULEu"))
  )
  expected_years_next_age <- matrix(
    0,
    nrow = 2,
    ncol = 2,
    dimnames = list(c("H", "U"), c("H", "U"))
  )
  
  for (age_index in rev(seq_len(n_ages))) {
    transient_matrix <- matrix(
      c(
        probability_wide$HH[age_index], probability_wide$HU[age_index],
        probability_wide$UH[age_index], probability_wide$UU[age_index]
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("H", "U"), c("H", "U"))
    )
    
    expected_years_this_age <- diag(2) + transient_matrix %*% expected_years_next_age
    expected_years[age_index, ] <- c(
      expected_years_this_age["H", "H"],
      expected_years_this_age["H", "U"],
      expected_years_this_age["U", "H"],
      expected_years_this_age["U", "U"]
    )
    expected_years_next_age <- expected_years_this_age
  }
  
  expectancy_by_age <- tibble(
    age = probability_wide$age,
    HLEh = expected_years[, "HLEh"],
    ULEh = expected_years[, "ULEh"],
    HLEu = expected_years[, "HLEu"],
    ULEu = expected_years[, "ULEu"]
  )
  
  result <- probability_wide %>%
    left_join(lxs, by = "age") %>%
    left_join(expectancy_by_age, by = "age") %>%
    mutate(
      HDt = HD * lh,
      HUt = HU * lh,
      UHt = UH * lu,
      UDt = UD * lu,
      dx = UDt + HDt,
      lx = lh + lu,
      ex = rev(cumsum(rev(lx))) / lx,
      edag = dx * ex,
      HLEdag_hu = HUt * (1 + lead(HLEh, default = 0) - lead(HLEu, default = 0)),
      HLEdag_hd = HDt * HLEh,
      HLEdag_uh = UHt * (-1 + lead(HLEu, default = 0) - lead(HLEh, default = 0)),
      HLEdag_ud = UDt * HLEu,
      ULEdag_hu = HUt * (-1 + lead(ULEh, default = 0) - lead(ULEu, default = 0)),
      ULEdag_hd = HDt * ULEh,
      ULEdag_uh = UHt * (1 + lead(ULEu, default = 0) - lead(ULEh, default = 0)),
      ULEdag_ud = UDt * ULEu,
      LEdag_hu = ULEdag_hu + HLEdag_hu,
      LEdag_hd = ULEdag_hd + HLEdag_hd,
      LEdag_uh = ULEdag_uh + HLEdag_uh,
      LEdag_ud = ULEdag_ud + HLEdag_ud
    ) 
  
  result
}

#' Calculate the 2D death distribution for one health, sex, year, and type.
#' Used by: R/09_manhattan.R.
calculate_group_dxh <- function(.x, .y,
                               expected_ages, required_transitions,
                               radix_age, last_age) {
  probability_wide <- .x %>%
    mutate(transition = paste0(as.character(from), as.character(to))) %>%
    select(age, transition, prob) %>%
    pivot_wider(names_from = transition, values_from = prob) %>%
    arrange(age)
  
  missing_transitions <- setdiff(required_transitions, names(probability_wide))
  if (length(missing_transitions) > 0) {
    stop(
      "A model group is missing transition(s): ",
      paste(missing_transitions, collapse = ", ")
    )
  }
  
  if (!identical(as.integer(probability_wide$age), expected_ages)) {
    stop(
      "Each model group must have one complete annual probability grid from age ",
      radix_age, " through ", last_age, "."
    )
  }
  
  if (anyNA(probability_wide[required_transitions])) {
    stop("A model group has an incomplete age-by-transition probability grid.")
  }
  
  initial_composition <- init_constant(probability_wide[1, ])
  
  # calc_dxh_fast() treats its final input row as the terminal age and forces
  # qH = qU = 1 there. This row closes the life table after the age-110 step.
  probability_wide <- probability_wide %>%
    add_row(
      age = last_age + 1L,
      HH = 0,
      HU = 0,
      HD = 1,
      UH = 0,
      UU = 0,
      UD = 1
    )
  
  calc_dxh_fast(
    p_tibble = probability_wide %>% select(age, all_of(required_transitions)),
    age_int = 1,
    init = initial_composition,
    keep_zero = FALSE
  )
}

