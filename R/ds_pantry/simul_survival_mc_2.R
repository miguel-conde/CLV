
# LIBRARIES & SOURCES -----------------------------------------------------


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

the_data <- list()

for (i in 1:N) {
  fst_1 <- rgeom(1, 1/24) + 1
  # print(fst_1)
  the_data[[i]] <- c(rep(0, fst_1-1), 1)
  names(the_data[[i]]) <- paste0("t_", 1:length(the_data[[i]]))
}


tbl_data <- bind_rows(the_data) %>% mutate_all(~ ifelse(is.na(.), 1, .))

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



mclistFit <- markovchainListFit(tbl_data, byrow = TRUE, laplacian = 0, name)

predict(mclistFit$estimate, c("0", "0"), n.ahead = 150, continue = FALSE)
predict(mclistFit$estimate, c("0", "0"), n.ahead = 100, continue = TRUE)

predict(mclistFit$estimate, c("0", "0"), n.ahead = 200, continue = FALSE)
predict(mclistFit$estimate, c("0", "0"), n.ahead = 200, continue = TRUE)
