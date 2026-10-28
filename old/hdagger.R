# packages
library(tidyverse)
library(collapse)
library(tidyfast)

# functions from Tim
source("functions_e_dagger.R")

# here are the transition probabilities 
# self reported health for 2013, single year, HMD adjustment
load("test_data_for_edagger.RData")


# filter only one year and one sex
ptibble <- final %>% 
  filter(time == 2013, sex == "female") %>% 
  unite(state, c("from", "to"), sep = "") %>%
  pivot_wider(names_from  = state,
              values_from = prob) %>% 
  dplyr::select(-c(1, 2)) |> 
  add_row(age=111,UD=1,HD=1,HH=0,HU=0,UH=0,UU=0)

# all good with sum constraints
ptibble %>% 
  dplyr::select(-1) %>% 
  rowSums()

# calculate initial constant
initial <- ptibble %>% 
  slice(1) %>%
  init_constant()

# calculate lh and lu
lh <- ptibble %>%
  p_tibble2lxs(init = initial, state = "H") %>%
  rename(lh = lxs)

lu <- ptibble %>% 
  p_tibble2lxs(init = initial, state = "U") %>%
  rename(lu = lxs)

# calculate the LE conditional on from state
LEi <- ptibble %>% 
    p_tibble2N() |>
    as.data.frame() |>
    rownames_to_column("state_age_to") |>
    pivot_longer(-state_age_to, names_to  = "state_age_from", values_to = "lxs") |>
    separate_wider_delim(
      cols  = state_age_from,
      delim = "_",
      names = c("state_from", "age_from")
    ) |>
    separate_wider_delim(
      cols  = state_age_to,
      delim = "_",
      names = c("state_to", "age_to")
    ) |>
    mutate(age_from = as.integer(age_from), 
           age_to   = as.integer(age_to)) |>
    filter(age_to >= age_from) |> 
    group_by(state_from, age_from, state_to) |>
    summarize(LEi = sum(lxs), .groups = "drop") |>
    mutate(expectancy = paste0(state_to, "LE", tolower(state_from))) |>
    select(!starts_with("state")) |>
    pivot_wider(names_from  = expectancy, values_from = LEi) |>
    rename(age = age_from)

# join and calculate daggers and lifetable for classic dagger
int <- 1 # ager interval

daggers <-
  lu %>% 
  full_join(lh, by = join_by(age)) %>%
  full_join(LEi, by = join_by(age)) %>%
  as_tibble() %>% 
  full_join(ptibble, by = join_by(age)) %>% 
  # calculate lt quantities
  mutate(lx = lh + lu,
         dx = lx - lead(lx, default = 0),
         ex = rev(cumsum(rev(lx))) / lx,
         edag = dx * ex) |> 
  filter(age <= 110) |> 
  # calculate transfers;
  # note UDt and HDt together sum to 1, which
  # means all transitions together > 1.
  mutate(HDt = HD * lh,
         HUt = HU * lh,
         UHt = UH * lu,
         UDt = UD * lu) |> 
  # t suffix stands for transfers
    mutate(
      # HLE
      HLEdag_hu = HUt * (int + lead(HLEh, n = int, default = 0) - lead(HLEu, n = int, default = 0)),
      HLEdag_hd = HDt * HLEh,
      HLEdag_uh = UHt * (-int + lead(HLEu, n = int,default = 0) - lead(HLEh, n = int,default = 0)),
      HLEdag_ud = UDt * HLEu,
      # ULE      
      ULEdag_hu = HUt * (-int + lead(ULEh, n = int,default = 0) - lead(ULEu, n = int,default = 0)),
      ULEdag_hd = HDt * ULEh,
      ULEdag_uh = UHt * (int + lead(ULEu, n = int,default = 0) - lead(ULEh, n = int,default = 0)),
      ULEdag_ud = UDt * ULEu,
      # LE
      LEdag_hu = ULEdag_hu + HLEdag_hu,
      LEdag_hd = ULEdag_hd + HLEdag_hd,
      LEdag_uh = ULEdag_uh + HLEdag_uh,
      LEdag_ud = ULEdag_ud + HLEdag_ud) |> 
    # w-suffix stands for weighted
    mutate(
      # HLE
      HLEdag_huw = HLEdag_hu / sum(HUt),
      HLEdag_hdw = HLEdag_hd / sum(HDt),
      HLEdag_uhw = HLEdag_uh / sum(UHt),
      HLEdag_udw = HLEdag_ud / sum(UDt),
      # ULE
      ULEdag_huw = ULEdag_hu / sum(HUt),
      ULEdag_hdw = ULEdag_hd / sum(HDt),
      ULEdag_uhw = ULEdag_uh / sum(UHt),
      ULEdag_udw = ULEdag_ud / sum(UDt),
      # LE (weights are the same for each row here)
      LEdag_huw = ULEdag_huw + HLEdag_huw,
      LEdag_hdw = ULEdag_hdw + HLEdag_hdw,
      LEdag_uhw = ULEdag_uhw + HLEdag_uhw,
      LEdag_udw = ULEdag_udw + HLEdag_udw) |> 
  # other quantities;
  # 'xd' means "by death" i.e. x is either h or u
  mutate(HLEdag_xdw = HLEh * HDt / sum(dx) +
                      HLEu * UDt / sum(dx),
         ULEdag_xdw = ULEh * HDt / sum(dx) +
                     ULEu * UDt / sum(dx),
         LEdag_xdw = HLEdag_xdw + ULEdag_xdw
         ) 
     
# survey results so far
 daggers |>
   select(age, contains("dag")) |> 
   summarize(across(-1,\(x) sum(x, na.rm = TRUE))) |> 
   pivot_longer(everything(), names_to = "variable", values_to = "value") |> 
   mutate(weighted = ifelse(grepl(variable,pattern="w"),"w","t"),
          variable = gsub(variable,pattern="w", replacement = "")) |> 
   pivot_wider(names_from = weighted, values_from = value) |> 
   View()
         
 daggers |>
   select(age, contains("dag")) |> 
   summarize(across(-1,\(x) sum(x, na.rm = TRUE))) |> 
   pivot_longer(everything(), names_to = "variable", values_to = "value") |> 
   filter(variable %in% c("LEdag_hu","LEdag_hd","LEdag_uh","LEdag_ud")) |> 
   pull(value) |> sum()
   
 daggers |>
   select(age, contains("dag")) |> 
   summarize(across(-1,\(x) sum(x, na.rm = TRUE))) |> 
   pivot_longer(everything(), names_to = "variable", values_to = "value") |> 
   filter(variable %in% c("LEdag_hd","LEdag_ud")) |> 
   pull(value) |> sum()

