
# initial constant
init_constant <- function(x) {
  x    <- unlist(x[c("HH", "UH", "HU", "UU")]) 
  u    <- matrix(x, nrow = 2, byrow = TRUE)
  v    <- eigen(u)$vectors[, 1]
  init <- v / sum(v)
  setNames(init, c("H", "U"))
}

# P2U
p_tibble2U <- function(p_tibble) {
  
  age <- p_tibble[["age"]] |>
    unique() |>
    sort()
  
  age <- c(age, age[length(age)] + 1)
  pre <- p_tibble |>
    dt_pivot_longer(-age, 
                    names_to  = "fromto", 
                    values_to = "p") |>
    fmutate(from = substr(fromto, 0, 1),
            to   = substr(fromto, 2, 2)) |>
    fselect(-fromto) |>
    fsubset(to != "D")
  
  states <- pre$from |>
    unique()
  
  pre |>
    group_by(from, to) |>
    nest() |>
    fmutate(data = map(data, ~ .x |>
                         pi_block_outer())) |> 
    dt_pivot_wider(names_from  = from, 
                   values_from = data) |> 
    unnest(cols = all_of(states),
           names_sep = "") |> 
    ungroup() |> 
    fmutate(to = paste(rep(states, each = length(age)), 
                       rep(age,    length(states)), # each
                       sep = "_")) |> 
    column_to_rownames("to") |> 
    as.matrix()
}

# P2N
p_tibble2N <- function(p_tibble, discount = FALSE) {
  
  U <- p_tibble2U(p_tibble)
  I <- diag(rep(1, nrow(U)))
  N <- solve(I - U) 
  if (discount) {
    N < N - I / 2
  }
  return(N)
}

# P2lx
p_tibble2lxs <- function(p_tibble, state = "H", init = c(H = 0.8, U = 0.2)) {
  
  age                 <- p_tibble[["age"]] |> min()
  N                   <- p_tibble2N(p_tibble)
  cols                <- grepl(colnames(N),
                               pattern =
                                 paste0("^._", age, "$"))
  rows                <- grepl(rownames(N), pattern = state)
  to_weight           <- N[rows, cols] 
  colnames(to_weight) <- substr(colnames(to_weight), 0, 1)  
  
  init <- as.matrix(init) |>
    as.data.frame() |>
    rownames_to_column("init_state") |>
    rename(init = V1)
  
  to_weight |>
    as.data.frame() |>
    rownames_to_column("age") |>
    as_tibble() |>
    fmutate(age = parse_number(age)) |>
    dt_pivot_longer(-age, 
                    names_to  = "init_state", 
                    values_to = "lxs") |>
    left_join(init, by = "init_state") |> 
    fgroup_by(age) |>
    fsummarize(lxs = sum(lxs * init)) |> 
    fungroup()
}

# helper
pi_block <- function(p, state_from, state_to, age) {
  
  state_fromi  <- state_from[1]
  state_toi    <- state_to[1]
  age          <- c(age, max(age) + 1)
  P            <- diag(p)
  P            <- cbind(rbind(0, P), 0)
  from_names   <- paste(state_fromi, age, sep = "_")
  to_names     <- paste(state_toi,   age, sep = "_")
  dimnames(P)  <- list(to_names, from_names)
  return(P)
}

# helper
pi_block_outer <- function(chunk) {
  
  pi_block(chunk[["p"]]|> as.double(),
           chunk[["from"]],
           chunk[["to"]],
           chunk[["age"]]) |>
    as_tibble()
}

# p_tibble2lhlu <- function(ptibble, init){
#   age <- ptibble$age
#   out <- tibble(age = age, 
#                 lh = c(init["H"],rep(0,length(age)-1)),
#                 lu = c(init["U"],rep(0,length(age)-1)))
#   for (i in 1:(nrow(out)-1)){
#     out$lh[i+1] <- ptibble$HH[i] * out$lh[i] + ptibble$UH[i] * out$lu[i]
#     out$lu[i+1] <- ptibble$UU[i] * out$lu[i] + ptibble$HU[i] * out$lh[i]
#   }
#   out
# }



