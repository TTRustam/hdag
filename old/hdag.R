
library(tidyverse)
source("https://raw.githubusercontent.com/timriffe/ms_dist/master/code/00_setup.R")
source("https://raw.githubusercontent.com/timriffe/ms_dist/master/code/01_functions.R")
adl <- read_csv("https://raw.githubusercontent.com/timriffe/ms_dist/master/data/share_all.csv",show_col_types = FALSE) |> 
  filter(measure == "ADL",
         age <=110) 
adl |> tail()
init <- init_constant(adl[1,])


lh <- 
  adl |>  
  select(-measure, -sex, -country, -version) |> 
  p_tibble2lxs(init = init, state="H") |> 
  rename(lh = lxs)
lu <- 
  adl |>  
  select(-measure, -sex, -country, -version) |> 
  p_tibble2lxs(init = init, state="U") |> 
  rename(lu = lxs)

LEi <- 
  adl |> 
  select(-measure, -sex, -country, -version) |> 
  p_tibble2N() |> 
  as.data.frame() |> 
  rownames_to_column("state_age_to") |> 
  pivot_longer(-state_age_to, 
               names_to = "state_age_from", 
               values_to = "lxs") |>   
  separate_wider_delim(cols = state_age_from, 
                       delim = "_",
                       names = c("state_from","age_from")) |> 
  separate_wider_delim(cols = state_age_to,
                       delim = "_",
                       names = c("state_to", "age_to")) |> 
  mutate(age_from = as.integer(age_from),
         age_to = as.integer(age_to)) |> 
  group_by(state_from, age_from, state_to) |> 
  summarize(LEi = sum(lxs),.groups = "drop") |> 
  mutate(expectancy = paste0(state_to,"LE",tolower(state_from))) |> 
  select(!starts_with("state")) |> 
  pivot_wider(names_from = expectancy, values_from = LEi) |> 
  rename(age = age_from)
                                                                    dags <-       
  adl |>   
  select(-HH, - UU) |>
  left_join(lh,by=join_by(age)) |> 
  left_join(lu,by=join_by(age)) |> 
  left_join(LEi, by = join_by(age)) |> 
  mutate(HLEdag_hu = HU * lh * (HLEh - HLEu),
         HLEdag_hd = HD * lh * HLEh,
         HLEdag_uh = UH * lu * (HLEu - HLEh), # negative
         HLEdag_ud = UD * lu * HLEu,
         
         ULEdag_hu = HU * lh * (ULEh - ULEu), # negative
         ULEdag_hd = HD * lh * ULEh,
         ULEdag_uh = UH * lu * (ULEu - ULEh),
         ULEdag_ud = UD * lu * ULEu,
         
         LEdag_hu = ULEdag_hu + HLEdag_hu,
         LEdag_hd = ULEdag_hd + HLEdag_hd,
         LEdag_uh = ULEdag_uh + HLEdag_uh,
         LEdag_ud = ULEdag_ud + HLEdag_ud) |> 
  select(age, contains("dag")) |> 
  pivot_longer(-age, 
               names_to = "dag_component", 
               values_to = "dag") |> 
  separate_wider_delim(dag_component, 
                       delim = "_", 
                       names = c("state_dag","transition")) 

                                                                                 dags |> 
  mutate(state_dag = factor(state_dag, 
                            levels = c("HLEdag","ULEdag","LEdag"))) |>
  ggplot(aes(x = age, y = dag, color = transition)) +
  geom_line() +
  facet_wrap(~state_dag) +
  theme_minimal()

dags |> 
  group_by(state_dag, transition) |> 
  summarise(dag = sum(dag)) |> 
  mutate(state_dag = factor(state_dag, 
         levels = c("HLEdag","ULEdag","LEdag"))) |>
  ggplot(aes(x=state_dag, fill = transition, y = dag)) +
  geom_col() +
  theme_minimal()

dags |> 
  group_by(state_dag) |> 
  summarize(dag = sum(dag))


# Here's a cheap lifetable edagger, low quality
n  = nrow(adl)
lx = lh$lh + lu$lu
dx = lh$lh * c(adl$HD,1) + lu$lu * c(adl$UD,1)
ex = rev(cumsum(rev(lx))) / lx
sum(dx * ex)

# compare ex:
prev = lh$lh / (lh$lh + lu$lu)
ex2 <- 
LEi |> 
  mutate(
    prev = prev,
    ex = prev * HLEh + prev * ULEh +
           (1-prev) * HLEu + (1-prev) * ULEu) |> 
  pull(ex)

LEi |> 
  ggplot(aes(x = age, y = HLEh)) +
  geom_line()

dags |> filter(state_dag == "LEdag") |> pull(dag) |> sum()
# this is smaller
dags |> 
  filter(transition %in% c("hd","ud"),
         state_dag == "LEdag") |> 
  summarize(sum(dag))

LEi |> 
  pivot_longer(-age, names_to = "type", values_to = "expectancy") |>   ggplot(aes(x = age, y = expectancy, color = type)) +
  geom_line() +
  theme_minimal()

# are these just sensitivities?
source("https://raw.githubusercontent.com/timriffe/ms_sensitivity/refs/heads/master/R/00_functions_classic.R")
source("https://raw.githubusercontent.com/timriffe/ms_sensitivity/refs/heads/master/R/00_sensitivity_functions.R")
spt <- 
p_tibble |> 
  pivot_longer(-age, names_to = "transition", values_to = "p") |> 
  s2t(init = init, expectancy = "h") |> 
  mutate(age = age + 50) |> 
  filter(transition !="init") 
spt |>   
ggplot(aes(x= age, y = effect, color = transition)) +
  geom_line() +
  theme_minimal()

s2t

