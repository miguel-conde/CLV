# https://danielweitzenfeld.github.io/passtheroc/blog/2015/01/19/s-b-g/

library(tidyverse)
library(rstan)
library(tidybayes)

df <- tribble(~cohorte, ~t, ~n,
              1, 0, 1000,
              1, 1, 840,
              1, 2, 618,
              1, 3, 570,
              1, 4, 519,
              2, 0, 1000,
              2, 1, 810,
              2, 2, 678,
              2, 3, 560,
              3, 0, 1000,
              3, 1, 900,
              3, 2, 790,
              4, 0, 1000,
              4, 1, 780)

df <- tribble(~cohorte, ~t, ~n,
              1, 0, 1000,
              1, 1, 869,
              1, 2, 743,
              1, 3, 653,
              1, 4, 593,
              1, 5, 551,
              1, 6, 517,
              1, 7, 491,
              1, 8, 468,
              1, 9, 445,
              1, 10, 427,
              1, 11, 409,
              1, 12, 394,
              2, 0, 1000,
              2, 1, 631,
              2, 2, 468,
              2, 3, 382,
              2, 4, 326,
              2, 5, 289,
              2, 6, 262,
              2, 7, 241,
              2, 8, 223,
              2, 9, 207,
              2, 10, 194,
              2, 11, 183,
              2, 12, 173)

model_data_all <- df %>% 
  group_by(cohorte) %>% 
  group_map(~ .x %>% arrange(t) %>% 
              mutate(lag_n = lag(n)) %>% 
              mutate(n_end = lag_n - n) %>% 
              mutate(r = n / lag_n ) %>% 
              drop_na() %>% 
              mutate(S = cumprod(r)) %>% 
              mutate(theta_emp = 1 - S^(1/t)) %>% 
              mutate(p_t_emp = (lag_n - n)/.$lag_n[1]), 
            .keep = TRUE) %>% 
  bind_rows()

# s-BG_0 ------------------------------------------------------------------

model_data_0 <- model_data_all %>% filter(cohorte == 1, t < 8)

STAN_FILE_0 <- here::here("stan", "s-BG_0.stan")
# STAN_FILE <- here::here("clv", "tst_stan.stan")

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

model_sbg_0 <- stan_model(file = STAN_FILE_0)

sbg_model_0 <- 
  sampling(model_sbg_0, 
           iter = 10000, chains = 1, 
           data = list(N = nrow(model_data_0),
                       n = model_data_0$n_end,
                       n_0 = model_data_0$lag_n[1],
                       mean_alpha = 1,
                       mean_beta = 1),
           verbose = TRUE)

sbg_model_0

model_data_0$r
extract(sbg_model_0)$sim_r %>% apply(2, median)

extract(sbg_model_0)$alpha %>% median()
extract(sbg_model_0)$beta %>% median()

extract(sbg_model_0)$alpha %>% hist()
extract(sbg_model_0)$beta %>% hist()

# s-BG --------------------------------------------------------------------

model_data <- model_data_all %>% filter(t < 8)

STAN_FILE <- here::here("clv", "s-BG.stan")
# STAN_FILE <- here::here("clv", "tst_stan.stan")

options(mc.cores = parallel:: detectCores())
rstan_options(auto_write = TRUE)

model_sbg <- stan_model(file = STAN_FILE)

sbg_model <- 
  sampling(model_sbg, 
           iter = 10000, chains = 1, 
           data = list(N = nrow(model_data),
                       N_g = model_data %>% count(cohorte) %>% pull(n),
                       G = 2,
                       g = model_data$cohorte,
                       n = model_data$n_end,
                       n_0 = c(1000, 1000)),
           verbose = TRUE)

sbg_model

check_hmc_diagnostics(sbg_model)
traceplot(sbg_model, c("alpha", "beta"))
stan_plot(sbg_model, c("alpha", "beta", "avg_alpha", "avg_beta"))
stan_hist(sbg_model, c("alpha", "beta", "avg_alpha", "avg_beta"))
stan_dens(sbg_model, c("alpha", "beta", "avg_alpha", "avg_beta"))

stan_diag(sbg_model, "sample")
stan_diag(sbg_model, "stepsize")
stan_diag(sbg_model, "treedepth")
stan_diag(sbg_model, "divergence")
stan_par(sbg_model, "alpha[1]")
stan_rhat(sbg_model, c("alpha", "beta", "avg_alpha", "avg_beta"))
stan_ess(sbg_model,  c("alpha", "beta", "avg_alpha", "avg_beta"))
stan_mcse(sbg_model, c("alpha", "beta", "avg_alpha", "avg_beta"))

model_data$r
extract(sbg_model)$sim_r %>% apply(2, median)

extract(sbg_model_0)$alpha %>% median()
extract(sbg_model_0)$beta %>% median()

extract(sbg_model_0)$alpha %>% hist()
extract(sbg_model_0)$beta %>% hist()


