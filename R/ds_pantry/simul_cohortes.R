library(tidyverse)
library(lubridate)

source(here::here("R", "src", "utils_cohorts.R"), encoding = "utf8")


# TEST UTILS --------------------------------------------------------------

# Cohortes por mes de cumpleaños ------------------------------------------


res_fun <- simul_births_deaths(N = 10000, 
                                vida_media = 60, 
                                f_gen_ini = "2003-01-01", f_gen_fin = "2022-12-31",
                                f_obs_ini = "2020-01-01", f_obs_fin = "2022-12-31", 
                                seed = 2023)

res_lifelines <- make_lifelines_month(res_fun, 2)

make_lifelines_month(res_fun, 12) %>% 
  select(id:yr_dth, ((ncol(.)-35):ncol(.))) %>% 
  align_lifelines() %>% 
  colSums(na.rm = TRUE) %>% 
  plot()

aligned_lifelines <- align_lifelines(res_lifelines)

make_emp_survivor_curve(aligned_lifelines)

make_lifelines_month(res_fun, mth = 1) %>% 
  align_lifelines() %>% 
  make_emp_survivor_curve() %>% 
  plot(type = "l")

make_lifelines_month(res_fun, mth = 3) %>% 
  align_lifelines() %>% 
  make_emp_survivor_curve() %>% 
  plot(type = "l")

make_lifelines_month(res_fun, mth = 12) %>% 
  align_lifelines() %>% 
  make_emp_survivor_curve() %>% 
  plot(type = "l")

make_lifelines_month(res_fun, mth = 9) %>% 
  align_lifelines() %>% 
  make_emp_survivor_curve()  %>% 
  make_emp_survivor_rate() %>% 
  plot(type = "l")

# Vida media empírica
make_lifelines_month(res_fun, mth = 3) %>% 
  align_lifelines() %>% 
  make_emp_survivor_curve(perc = TRUE) %>% sum()

# Cohortes por mes de nacimiento ------------------------------------------

# Cohorte de los nacidos en los meses de febrero en los años de observación
make_lifelines_month(res_fun, 2) %>% 
  filter(yr_bth %in% 2020:2022) %>% 
  align_lifelines()%>% 
  make_emp_survivor_curve() %>% 
  plot(type = "l")

make_lifelines_month(res_fun, 2) %>% 
  filter(yr_bth %in% 2020:2022) %>% 
  align_lifelines()%>% 
  make_emp_survivor_curve() %>% 
  make_emp_survivor_rate() %>% 
  plot(type = "l")

# ORIGINAL ----------------------------------------------------------------



# Simulate cohorts

N <- 100000 # Nº de individuos

# Duración de la vida - Distr. geométrica, vida media = 60
p = 1/60

# Tiempo discreto

aux_dates <- tibble(t = 1:240,
                    date = seq.Date(as_date("2003-01-01"),
                                    as_date("2022-12-31"),
                                    by = "month")) %>%
  mutate(mth = lubridate::month(date),
         yr = lubridate::year(date))



# Los nacimientos se distribuyen uniformemente
set.seed(2023)
probe <- tibble(t_bth = sample(1:240, size = N, replace = TRUE),
                dur_vida = rgeom(N, p),
                t_dth = t_bth + dur_vida) %>%
  arrange(t_bth) %>%
  left_join(aux_dates, by = c("t_bth" = "t")) %>%
  mutate(date_dth = date + months(dur_vida)) %>%
  mutate(mth_dth = lubridate::month(date_dth, label = TRUE, abbr = FALSE),
         yr_dth = lubridate::year(date_dth)) %>%
  select(t_bth, dur_vida, t_dth,
         date_bth = date, mth_bth = mth, yr_bth = yr,
         date_dth, mth_dth, yr_dth) %>% 
  # Cojo los activos en el periodo 2020-2022
  filter(yr_dth > 2019) %>%
  mutate(date_dth = ifelse(date_dth >= as_date("2023-01-01"), NA, date_dth) %>% 
           as_date()) %>% 
  mutate(mth_dth  = ifelse(date_dth >= as_date("2023-01-01"), NA, mth_dth)) %>%
  mutate(yr_dth   = ifelse(date_dth >= as_date("2023-01-01"), NA, yr_dth)) %>%
# mutate(dur_vida   = ifelse(date_dth >= as_date("2023-01-01"), NA, dur_vida))
  mutate(id = row_number(), .before = 1)

probe %>%
  group_by(mth_bth) %>%
  summarise(n = n(), avg_dur_vida = mean(dur_vida, na.rm = TRUE))

probe %>%
  group_by(mth_dth) %>%
  summarise(n = n(), avg_dur_vida = mean(dur_vida, na.rm = TRUE))


kk <- purrr::map2(probe %>%
                    filter(mth_bth == 1) %>%
                    pull(date_bth),
                  probe %>%
                    filter(mth_bth == 1) %>%
                    pull(date_dth),
                  function(x,y)  {
                    if (is.na(y)) y <- as_date("2022-12-31")
                    the_names <- seq.Date(x, y, by = "month")
                    out <- rep(1, length(the_names))
                    names(out) <- the_names
                    return(out)
                  })

res <- bind_cols(probe %>% filter(mth_bth == 1)) %>% 
  bind_cols(kk %>% 
              bind_rows())

res %>%
  select(-(id:yr_dth)) %>% 
  summarise_all(~ sum(., na.rm = TRUE))  %>%
  gather(date, obs) %>%
  mutate(date = as_date(date)) %>% 
  filter(date > as_date("2019-12-31")) %>% 
  plot(type = "l", ylim = c(0, 300))

res

res2 <- res %>% 
  select(-(id:yr_dth)) %>% 
  rowSums(na.rm = TRUE) %>% 
  map(~ {
    out <- rep(1, .x)
    names(out) <- paste0("mth_", 1:length(out))
    return(out)
    })  %>% 
  bind_rows() %>% 
  colSums(na.rm = TRUE)

res %>% 
  select(-(id:yr_dth)) %>% 
  rowSums(na.rm = TRUE) %>% 
  table %>% 
  as_tibble() %>% 
  rename(age_mths = 1) %>% 
  mutate(age_mths = as.integer(age_mths)) %>% 
  arrange(desc(age_mths)) %>% 
  plot()

res %>% 
  select(-(id:yr_dth)) %>% 
  rowSums(na.rm = TRUE) %>% 
  hist()

# Función supervivencia
(res2/lag(res2))[-1] %>% cumprod() %>% plot
curve((1-p)^x, 1, 240, add = TRUE)

bayes_data <- (-diff(res2))
sbg_model_0 <- 
  sampling(model_sbg_0, 
           iter = 1000, chains = 1, 
           data = list(N = length(bayes_data),
                       n = bayes_data,
                       n_0 = res2[1],
                       mean_alpha = 1,
                       mean_beta = 1),
           verbose = TRUE)

sbg_model_0

(res2/lag(res2))[-1]
extract(sbg_model_0)$sim_r %>% apply(2, median)

plot((res2/lag(res2))[-1],
     extract(sbg_model_0)$sim_r %>% apply(2, median))

extract(sbg_model_0)$alpha %>% median()
extract(sbg_model_0)$beta %>% median()

extract(sbg_model_0)$alpha %>% hist()
extract(sbg_model_0)$beta %>% hist()


library(BTYD)
BTYD::bgbb.EstimateParameters()