
# LIBRARIES & SOURCES -----------------------------------------------------
library(tidyverse)
library(markovchain)

N <- 1000


# FUNCTIONS ---------------------------------------------------------------

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

bs_fvals <- function(tbl_data, initialstate = c(1,0), n = 36,
                     n_bs = 100, mode = c("freqs", "probs"), seed = NULL) {
  
  mode <- match.arg(mode, choices = c("freqs", "probs"))
  probs <- ifelse(mode == "probs", TRUE, FALSE)
  
  n_r <- nrow(tbl_data)
  
  lst_muertos <- list()
  states_0_1 <- c("0", "1")
  byRow <- TRUE
  
  set.seed(seed)
  for (i in 1:n_bs) {
    muertos_i <- tbl_data[sample(1:n_r, n_r, replace = TRUE), ] %>% 
      colSums()
    
    tbl_survival <- bind_rows(tibble(vivos = n_r, muertos = 0,),
                              tibble(vivos = n_r - muertos_i,
                                     muertos = muertos_i)) %>% 
      mutate(i = 0:(nrow(.)-1), .before = 1)
    
    lst_transitions_matrix <- list()
    
    for (j in 2:nrow(tbl_survival)) {
      
      # Survivor function
      surv_fun <- as.numeric(tbl_survival[j, "vivos"]) / N
      matrix_0_1 <- matrix(data = c(surv_fun, 1 - surv_fun,  0,  1) , 
                           byrow = byRow, 
                           nrow = 2,
                           dimnames = list(states_0_1, states_0_1))
      
      lst_transitions_matrix[[j-1]] <- new(Class            = "markovchain", 
                                           states           = states_0_1, 
                                           byrow            = byRow,
                                           transitionMatrix = matrix_0_1, 
                                           name             = "De 0 a 1")
      
    }
    
    mc_propensity <- new("markovchainList",
                         markovchains = lst_transitions_matrix,
                         name = "Propensity")
    
    lst_muertos[[i]] <- fvals(mchain = mc_propensity, initialstate = initialstate, n = n)
  }
  
  aux <- lst_muertos[[1]]
  for (i in 2:length(lst_muertos)) {
    aux <- abind::abind(aux, lst_muertos[[i]], along = 3)
  }
  
  out <- list(avgs = aux %>% apply(c(1,2), mean),
              sds  = aux %>% apply(c(1,2), sd))
  
  return(out)
    
}

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





# DATA --------------------------------------------------------------------

make_time_to_event_dataset <- function(N_t, N, p, states = c("A", "D"),
                                       seed = NULL, format = c("list", "data.frame")) {
  
  format <- match.arg(format, c("list", "data.frame"))
  
  the_data <- list()
  
  set.seed(seed)
  for (i in 1:N) {
    fst_1 <- rgeom(1, p) + 1
    the_data[[i]] <- c(rep(states[1], fst_1-1), states[2])[1:N_t]
    the_data[[i]][is.na(the_data[[i]])] <- states[2]
    names(the_data[[i]]) <- paste0("t_", 1:length(the_data[[i]]))
  }
  
  if (format == "list") return(the_data)
  
  tbl_data <- bind_rows(the_data) %>% mutate_all(~ ifelse(is.na(.), "D", .))
  
  return(tbl_data)
}

make_time_to_event_dataset(N_t = 36, N = N, p = 1/24)
make_time_to_event_dataset(N_t = 36, N = N, p = 1/24, format = "data.frame")
make_time_to_event_dataset(N_t = 36, N = N, p = 1/24, 
                           states = c("Alive", "Death"), format = "data.frame")

the_data <- list()
N_obs <- 36

for (i in 1:N) {
  fst_1 <- rgeom(1, 1/24) + 1
  # print(fst_1)
  the_data[[i]] <- c(rep("A", fst_1-1), "D")[1:N_obs]
  the_data[[i]][is.na(the_data[[i]])] <- "D"
  names(the_data[[i]]) <- paste0("t_", 1:length(the_data[[i]]))
}


tbl_data <- bind_rows(the_data) %>% mutate_all(~ ifelse(is.na(.), "D", .))

# CHICHA ------------------------------------------------------------------

vivos_muertos <- bs_tbl_data(tbl_data)

tbl_survival <- bind_rows(tibble(vivos = N, muertos = 0,),
                          tibble(vivos = vivos_muertos$vivos,
                                 muertos = vivos_muertos$muertos)) %>% 
  mutate(i = 0:(nrow(.)-1), .before = 1)



res <- bs_fvals(tbl_data, initialstate = c(1,0), n = 36,
                n_bs = 100, mode = "probs", seed = NULL)

plot(res$avgs[, "0"], type = "l")
lines((N - colSums(tbl_data))[1:36] / N, col = "red")

plot((N - colSums(tbl_data))[1:36], res$avgs[, "0"])


# CHICHA 2 ----------------------------------------------------------------

matrix_V_M <- createSequenceMatrix(stringchar = the_data %>% lapply(as.character), 
                                   toRowProbs = TRUE)
states_V_M <- c("Muertos", "Vivos")
dimnames(matrix_V_M) <- list(c(), states_V_M)

mcFitMLE <- markovchainFit(data = the_data %>% lapply(as.character))
mcFitBSP <- markovchainFit(data = tbl_data, 
                           method = "bootstrap", nboot = 5, name = "Bootstrap Mc")



mclistFit <- markovchainListFit(tbl_data, byrow = TRUE, laplacian = 0, name = "List MC")

predict(mclistFit$estimate, c("A", "A"), n.ahead = 150, continue = FALSE)
predict(mclistFit$estimate, c("D", "D"), n.ahead = 100, continue = TRUE)

predict(mclistFit$estimate, c("0", "0"), n.ahead = 200, continue = FALSE)
predict(mclistFit$estimate, c("0", "0"), n.ahead = 200, continue = TRUE)


# PREDIC ------------------------------------------------------------------

my_mclist_predict <- function(mclistFit_estimate, seq_prev, n_ahead) {
  
  fst_pred_state <- length(seq_prev)
  lst_pred_state <- min(n_ahead+fst_pred_state-1, 
                        length(mclistFit_estimate@markovchains))
  
  act_state <- vector(mode = "list", length = lst_pred_state - fst_pred_state + 2)
  act_state[[1]] = as.numeric(mclistFit_estimate[[1]]@states == tail(seq_prev,1))
  names(act_state[[1]]) <- mclistFit_estimate[[1]]@states
  
  for (i in fst_pred_state:lst_pred_state) {
    act_state[[i]] <- act_state[[i-1]] * mclistFit_estimate@markovchains[[i]]
  }
  
  out <- act_state[-1] %>% 
    lapply(as.data.frame) %>% 
    bind_rows() %>% 
    mutate(t = fst_pred_state:lst_pred_state, .before = 1) %>% 
    rowwise() %>% 
    mutate(pred = mclistFit_estimate[[1]]@states[which(c_across(2:ncol(.)) == max(c_across(2:ncol(.))))[1]]) %>% 
    ungroup()
  
  return(out)
}


my_mclist_predict(mclistFit$estimate, seq_prev = c("A", "A"), n_ahead = 18)
predict(mclistFit$estimate, newdata = c("A", "A"), n.ahead = 18)

rmarkovchain(20, mcCCRC, what = "list")

predict(mcCCRC, newdata = c("H", "H"), n.ahead = 5)
my_mclist_predict(mcCCRC, seq_prev = c("H", "H"), n_ahead = 5)


kk <- tibble(A = sapply(tbl_data, function(x) sum( x=="A")), 
             B = sapply(tbl_data, function(x) sum(x == "D"))) %>%
  mutate(rr = c(.$A[1] / (.$A[1]+.$B[1]), (A/lag(A))[-1]), 
         churn_rate = 1 - rr, 
         survivor_fun = cumprod(rr))
kk

kk$rr %>% plot(type = "l")
kk$rr %>% density() %>% plot()

kk$churn_rate %>% plot(type = "l")
kk$churn_rate %>% density() %>% plot()

kk$survivor_fun %>% plot(type = "l")

kk$survivor_fun %>% sum()

