library(tidyverse)
library(markovchain)

N <- 1000


# SURVIVAL ----------------------------------------------------------------

the_data <- list()

the_lengths <- rgeom(N, 1/24) + 1

for (i in 1:N) {
  fst_1 <- sample(1:100, 1)
  the_data[[i]] <- c(rep(0, fst_1-1), rep(1, the_lengths[i])) 
  names(the_data[[i]]) <- paste0("t_", 1:length(the_data[[i]]))
}

tbl_data <- bind_rows(the_data) %>% 
  mutate_all(~ ifelse(is.na(.), 0, .)) %>% 
  t()

colnames(tbl_data) <- paste0("C_", 1:N)
tbl_data <- tbl_data %>% t() %>% as_tibble()

de_0_a_1 <- tbl_data %>% t() %>% as_tibble() %>% 
  mutate_all(~ (c(NA, diff(.)) == -1)) %>% 
  sum(na.rm = TRUE)
de_1_a_0 <- tbl_data %>% t() %>% as_tibble() %>% 
  mutate_all(~ (c(NA, diff(.)) ==  1)) %>% 
  sum(na.rm = TRUE)
de_1_a_1 <- tbl_data %>% t() %>% as_tibble() %>% 
 summarise_all(~ sum(.[. == lag(.)] == 1, na.rm = TRUE)) %>% 
  sum()
de_0_a_0 <- tbl_data %>% t() %>% as_tibble() %>% 
  summarise_all(~ sum(.[. == lag(.)] == 0, na.rm = TRUE)) %>% 
  sum()

de_all <- de_0_a_1 + de_1_a_0 + de_1_a_1 + de_0_a_0


states_0_1 <- c("Non Contract", "Contract")
byRow <- TRUE
matrix_0_1 <- matrix(data = c(de_0_a_0 / (de_0_a_0 + de_0_a_1), de_0_a_1 / (de_0_a_0 + de_0_a_1),
                              de_1_a_0 / (de_1_a_0 + de_1_a_1), de_1_a_1 / (de_1_a_0 + de_1_a_1)), 
                     byrow = byRow, 
                     nrow = 2,
                     dimnames = list(states_0_1, states_0_1))

matrix_0_1 <- createSequenceMatrix(stringchar = the_data %>% lapply(as.character), 
                                   toRowProbs = TRUE)
dimnames(matrix_0_1) <- list(states_0_1, states_0_1)

mc_de_0_a_1 <- new(Class            = "markovchain", 
                   states           = states_0_1, 
                   byrow            = byRow,
                   transitionMatrix = matrix_0_1, 
                   name             = "De 0 a 1")

fvals <- function(mchain, initialstate, n) {
  
  out <- data.frame()
  names(initialstate) <- names(mchain)
  for (i in 0:n)
  {
    iteration <- initialstate*mchain^(i)
    out <- rbind(out, iteration)
  }
  out <- cbind(out, i = seq(0, n))
  out <- out[, c(3, 1:2)]
  
  return(out)
}

initialState <- c(N,0)
initialState * mc_de_0_a_1^10

fvals(mchain = mc_de_0_a_1, initialstate = initialState, n = 4)

rmarkovchain(100, mc_de_0_a_1, t = "Non Contract")
markovchainFit(matrix_0_1)

# PROPENSIÓN --------------------------------------------------------------

the_data <- list()

for (i in 1:N) {
  fst_1 <- rgeom(1, 1/24) + 1
  # print(fst_1)
  the_data[[i]] <- c(rep(0, fst_1-1), 1)
  names(the_data[[i]]) <- paste0("t_", 1:length(the_data[[i]]))
}


tbl_data <- bind_rows(the_data) %>% mutate_all(~ ifelse(is.na(.), 1, .))

bs_tbl_data <- function(tbl_data, 
                        n = 100, mode = c("freqs", "probs"), seed = NULL) {
  
  mode <- match.arg(mode, choices = c("freqs", "probs"))
  probs <- ifelse(mode == "probs", TRUE, FALSE)
  
  n_r <- nrow(tbl_data)
  
  aux <- list()
  
  set.seed(seed)
  for (i in 1:n) {
    aux[[i]] <- tbl_data[sample(1:n_r, n_r, replace = TRUE), ] %>% 
      colSums()
  }
  aux <- aux %>% bind_rows()
  
  out <- tibble(vivos = n_r - colMeans(aux),
                muertos = n_r - vivos) %>% 
    mutate(sd = apply(aux, 2, function(x) sd(x)))
  
  if (probs == TRUE) {
    out <- out %>% 
      mutate(vivos = vivos / n_r,
             muertos = muertos / n_r,
             sd = sd / n_r)
  }
  
  return(out)
}

bs_fvals <- function(mchain = mc_propensity, initialstate = c(1,0), n = 36) {
  
  bs_tbl_data(tbl_data, 
                          n = 100, mode = c("freqs", "probs"), seed = NULL)
  
}

vivos_muertos <- bs_tbl_data(tbl_data)

tbl_survival <- bind_rows(tibble(vivos = N, muertos = 0,),
                          tibble(vivos = vivos_muertos$vivos,
                                 muertos = vivos_muertos$muertos)) %>% 
  mutate(i = 0:(nrow(.)-1), .before = 1)


states_0_1 <- c("0", "1")
byRow <- TRUE
lst_transitions_matrix <- list()

for (i in 2:nrow(tbl_survival)) {
  
  # Survivor function
  surv_fun <- as.numeric(tbl_survival[i, "vivos"]) / N
  matrix_0_1 <- matrix(data = c(surv_fun, 1 - surv_fun,  0,  1) , 
                       byrow = byRow, 
                       nrow = 2,
                       dimnames = list(states_0_1, states_0_1))
  
  lst_transitions_matrix[[i-1]] <- new(Class            = "markovchain", 
                                       states           = states_0_1, 
                                       byrow            = byRow,
                                       transitionMatrix = matrix_0_1, 
                                       name             = "De 0 a 1")
  
}

mc_propensity <- new("markovchainList",
                     markovchains = lst_transitions_matrix,
                     name = "Propensity")

fvals <- function(mchain, initialstate, n) {
  
  out <- data.frame()
  names(initialstate) <- names(mchain)
  for (i in 1:n)
  {
    iteration <- initialstate*mchain[[i]]
    out <- rbind(out, iteration)
  }
  out <- cbind(out, i = seq(1, n))
  out <- out[, c(3, 1:2)]
  
  return(out)
}

initialState <- c(N, 0)

fvals(mchain = mc_propensity, initialstate = initialState, n = 4)

evol <- fvals(mchain = mc_propensity, initialstate = c(1,0), n = 36)

survived <- evol[, 2]
survival_rate <- survived / lag(survived)

survival_rate[-1] %>% plot(type = "l")

plot(survived)
plot(tbl_survival$vivos[2:37], survived)

steadyStates(mc_propensity[[1]])
plot(mc_propensity[[1]])

# We conclude that the probability for new contract  to be the third one, 
# given that the current state is non-contract, is given by:
firstPassagePdF <- firstPassage(object = mc_propensity[[1]], state = "0", n = 100)
firstPassagePdF[3, "1"]

firstPassagePdF[, "1"] %>% plot(type = "l")

meanFirstPassageTime(mc_propensity[[1]])

firstPassagePdF.long <- firstPassage(object = mc_propensity[[1]], state = "0", n = 100)
sum(firstPassagePdF.long[,"1"] * 1:100)


### Fit


# REVOLVING ---------------------------------------------------------------

library(lubridate)
source(here::here("R", "src", "utils_cohorts.R"), encoding = "utf8")

simul_feb <- simul_N("2020-02-01", "2022-12-31", N = 1000)

tbl_data <- simul_feb %>% spread(date, sales) %>% select(-id)

de_0_a_1 <- tbl_data %>% 
  mutate_all(~ (c(NA, diff(.)) == -1)) %>% 
  sum(na.rm = TRUE)
de_1_a_0 <- tbl_data %>% 
  mutate_all(~ (c(NA, diff(.)) ==  1)) %>% 
  sum(na.rm = TRUE)
de_1_a_1 <- tbl_data %>% 
  summarise_all(~ sum(.[. == lag(.)] == 1, na.rm = TRUE)) %>% 
  sum()
de_0_a_0 <- tbl_data %>% 
  summarise_all(~ sum(.[. == lag(.)] == 0, na.rm = TRUE)) %>% 
  sum()

states_0_1 <- c("0", "1")
byRow <- TRUE
matrix_0_1 <- matrix(data = c(de_0_a_0 / (de_0_a_0 + de_0_a_1), de_0_a_1 / (de_0_a_0 + de_0_a_1),
                              de_1_a_0 / (de_1_a_0 + de_1_a_1), de_1_a_1 / (de_1_a_0 + de_1_a_1)), 
                     byrow = byRow, 
                     nrow = 2,
                     dimnames = list(states_0_1, states_0_1))

mc_de_0_a_1 <- new(Class            = "markovchain", 
                   states           = states_0_1, 
                   byrow            = byRow,
                   transitionMatrix = matrix_0_1, 
                   name             = "De 0 a 1")

initialState <- c(N,0)
initialState * mc_de_0_a_1^10
fvals(mchain = mc_de_0_a_1, initialstate = initialState, n = 4)

fvals(mchain = mc_de_0_a_1, initialstate = initialState, n = 36)
