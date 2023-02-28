
simul_births_deaths <- function(N = 10000, 
                                vida_media = 60, 
                                f_gen_ini = "2003-01-01", f_gen_fin = "2022-12-31",
                                f_obs_ini = "2020-01-01", f_obs_fin = "2022-12-31",
                                seed = NULL) {
  # N: Nº de individuos
  # vida_media: inverso de la probabilidad de muerte
  # f_gen_ini, f_gen_fin: fechas iniciales y finales de generacion
  
  require(lubridate)
  
  f_gen_ini <- as_date(f_gen_ini)
  f_gen_fin <- as_date(f_gen_fin)
  f_obs_ini <- as_date(f_obs_ini)
  f_obs_fin <- as_date(f_obs_fin)
  
  p <- 1/vida_media
  
  seq_dates <- seq.Date(f_gen_ini, f_gen_fin, by = "month")
  
  # T: periodo de generación
  T <- length(seq_dates)
  
  # Tiempo discreto
  aux_dates <- tibble(t = 1:T, date = seq_dates) %>%
    mutate(mth = lubridate::month(date),
           yr = lubridate::year(date))
  
  set.seed(seed)
  # Los nacimientos se distribuyen uniformemente
  out <- tibble(t_bth = sample(1:240, size = N, replace = TRUE),
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
    filter(date_dth >= f_obs_ini) %>%
    mutate(date_dth = ifelse(date_dth > f_obs_fin, NA, date_dth) %>% 
             as_date()) %>% 
    mutate(mth_dth  = ifelse(date_dth > f_obs_fin, NA, mth_dth)) %>%
    mutate(yr_dth   = ifelse(date_dth > f_obs_fin, NA, yr_dth)) %>%
    mutate(id = row_number(), .before = 1)
  
  return(out)
  
}

make_lifelines_month <- function(probe, mth = 1) {
  
  max_date <- max(probe %>%
                    filter(mth_bth == mth) %>%
                    pull(date_dth), na.rm = TRUE)
  
  aux <- purrr::map2(probe %>%
                       filter(mth_bth == mth) %>%
                       pull(date_bth),
                     probe %>%
                       filter(mth_bth == mth) %>%
                       pull(date_dth),
                     function(x,y)  {
                       if (is.na(y)) y <- max_date
                       the_names <- seq.Date(x, y, by = "month")
                       out <- rep(1, length(the_names))
                       names(out) <- the_names
                       return(out)
                     })
  
  out <- probe %>% filter(mth_bth == mth) %>% 
    bind_cols(aux %>% 
                bind_rows())
  
  return(out)
}

# Cohortes por mes de nacimiento ------------------------------------------

make_birth_cohort <- function(N = 10000, 
                              vida_media = 60, 
                              f_gen_ini = "2003-01-01", f_gen_fin = "2022-12-31",
                              f_obs_ini = "2020-01-01", f_obs_fin = "2022-12-31",
                              seed = NULL) {
  
  
}

# Cohortes por mes de cumpleaños ------------------------------------------



align_lifelines <- function(probe) {
  
  out <- probe %>% 
    select(-(id:yr_dth)) %>% 
    rowSums(na.rm = TRUE) %>% 
    map(~ {
      out <- rep(1, .x)
      names(out) <- paste0("mth_", 1:length(out))
      return(out)
    })  %>% 
    bind_rows() # %>% 
  # colSums(na.rm = TRUE)
  
  return(out)
  
}

make_emp_survivor_curve <- function(probe, perc = FALSE) {
  
  out <- probe %>% 
    colSums(na.rm = TRUE)
  
  if (perc == TRUE) out <- out / out[1]
  
  return(out)
}

make_emp_survivor_rate <- function(probe) {
  
  out <- (probe / lag(probe))[-1]
  
  return(out)
}



# Simulación revolving ----------------------------------------------------

# Para cada cliente simulamos:
#
# - Duración de la primera activación - Geométrica(p_0) - Media 12
# - Número de reactivaciones - Poisson(lambda) - Media 5
# - Para cada reactivación i:
#   - Duración del intervalo inicial inactivo - Geométrica(p_i_i) - Media 3
#   - Duración del intervalo activo - Geométrica(p_a_i) - Media 6
#
#
# - p_0, p_i_i y p_a_i - beta(alpha, beta)

simul_1 <- function(p_0 = 1/12, lambda = 5, p_i_i = 1/3, p_a_i = 1/6, 
                    seed = NULL) {
  
  set.seed(seed)
  
  # Duración de la primera activación - Geométrica(p_0) - Media 12
  a_p_0 <- rbeta(1, shape1 = 10, shape2 = 10/p_0 - 10)
  # duracion_activacion_1 <- rgeom(1, prob = a_p_0)
  duracion_activacion_1 <- max(rnbinom(1, size = 1, prob = a_p_0), 1)
  
  out <- rep(1, duracion_activacion_1)
  
  # Número de reactivaciones - Poisson(lambda) - Media 5
  a_lambda <- rexp(1, 1/lambda)
  n_reactivaciones <- rpois(1, a_lambda)
  
  for (r in 1:n_reactivaciones) {
    # Duración del intervalo inicial inactivo - Geométrica(p_i_i) - Media 3
    a_p_i_r <- rbeta(1, shape1 = 10, shape2 = 10/p_i_i - 10)
    duracion_inactividad_r <- rgeom(1, prob = a_p_i_r)
    
    out <- c(out, rep(0, duracion_inactividad_r))
    
    # Duración del intervalo activo - Geométrica(p_a_i) - Media 6
    a_p_a_r <- rbeta(1, shape1 = 10, shape2 = 10/p_a_i - 10)
    # duracion_reactividad_r <- max(rgeom(1, prob = a_p_a_r), 1)
    duracion_reactividad_r <- max(rnbinom(1, size = 1, prob = a_p_a_r), 1)
    
    out <- c(out, rep(1, duracion_reactividad_r))
  }
  
  return(out)
}

simul_N <- function(fecha_ini, fecha_fin, 
                    N = 100, p_0 = 1/12, lambda = 5, p_i_i = 1/3, p_a_i = 1/6, 
                    seed = NULL) {
  
  fecha_ini <- as_date(fecha_ini)
  fecha_fin <- as_date(fecha_fin)
  
  out <- lapply(1:N, function(x) {
    aux <- simul_1()
    the_names <- seq.Date(fecha_ini, fecha_ini + months(length(aux)-1), by = "month")
    names(aux) <- the_names
    
    return(aux)
  }) %>% 
    bind_rows() %>% 
    select(as.character(seq.Date(fecha_ini, fecha_fin, by = "month"))) %>% 
    mutate(id = row_number(), .before = 1) %>% 
    gather(date, sales, -id) %>% 
    arrange(id) %>% 
    mutate(date = as_date(date))
  
  return(out)
}

build_rf_table <- function(probe) {
  
  aux <- probe %>% spread(date, sales) %>% select(-id) 
  
  out <- tibble(f = aux %>% rowSums(na.rm = TRUE),
                r = aux %>% apply(1, function(x) max(which(x == 1)))) %>% 
    group_by(f, r) %>% 
    summarise(n = n(), .groups = "drop") %>% 
    arrange(desc(f), desc(r)) %>% 
    mutate(T = ncol(aux))
  
  return(out)
}


# Bug BYTD ----------------------------------------------------------------

my_dc.MakeRFmatrixCal <- 
  function (frequencies, periods.of.final.purchases, num.of.purchase.periods, 
            holdout.frequencies = NULL) 
{
  if (!is.numeric(periods.of.final.purchases)) {
    stop("periods.of.final.purchases must be numeric")
  }
  if (length(periods.of.final.purchases) != length(frequencies)) {
    stop(paste("number of customers in frequencies is not equal", 
               "to the last purchase period vector"))
  }
  rf.mx.skeleton <- dc.MakeRFmatrixSkeleton(num.of.purchase.periods)
  if (is.null(holdout.frequencies)) {
    RF.matrix <- cbind(rf.mx.skeleton, num.of.purchase.periods, 
                       0)
    colnames(RF.matrix) <- c("x", "t.x", "n.cal", "custs")
  }
  else {
    RF.matrix <- cbind(rf.mx.skeleton, num.of.purchase.periods, 
                       0, 0)
    colnames(RF.matrix) <- c("x", "t.x", "n.cal", "custs", 
                             "x.star")
  }
  rf.n.custs <- cbind(frequencies, periods.of.final.purchases, 
                      holdout.frequencies)
  zeroes.rf.subset <- which(rf.n.custs[, 1] == 0)
  RF.matrix[1, 4] <- length(zeroes.rf.subset)
  if (!is.null(holdout.frequencies)) {
    RF.matrix[1, 5] <- sum(holdout.frequencies[zeroes.rf.subset])
  }
  # rf.n.custs <- rf.n.custs[-zeroes.rf.subset, ]
  if (length(zeroes.rf.subset) > 0) rf.n.custs <- rf.n.custs[-zeroes.rf.subset, ]
  rf.n.custs <- rf.n.custs[order(rf.n.custs[, 1], rf.n.custs[, 
                                                             2]), ]
  current.pair <- c(rf.n.custs[1, 1], rf.n.custs[1, 2])
  same.item.in.a.row.counter <- 1
  if (!is.null(holdout.frequencies)) {
    x.star.total <- rf.n.custs[1, 3]
  }
  num.count.points <- nrow(rf.n.custs)
  for (ii in 2:num.count.points) {
    last.pair <- current.pair
    current.pair <- c(rf.n.custs[ii, 1], rf.n.custs[ii, 2])
    if (identical(last.pair, current.pair)) {
      same.item.in.a.row.counter <- same.item.in.a.row.counter + 
        1
      if (!is.null(holdout.frequencies)) {
        x.star.total <- x.star.total + rf.n.custs[ii, 
                                                  3]
      }
    }
    else {
      x <- last.pair[1]
      t.x <- last.pair[2]
      corresponding.rf.index <- (x - 1) + 1 + t.x * (t.x - 
                                                       1)/2 + 1
      RF.matrix[corresponding.rf.index, 4] <- same.item.in.a.row.counter
      same.item.in.a.row.counter <- 1
      if (!is.null(holdout.frequencies)) {
        RF.matrix[corresponding.rf.index, 5] <- x.star.total
        x.star.total <- rf.n.custs[ii, 3]
      }
    }
    if (ii == num.count.points) {
      x <- current.pair[1]
      t.x <- current.pair[2]
      corresponding.rf.index <- (x - 1) + 1 + t.x * (t.x - 
                                                       1)/2 + 1
      RF.matrix[corresponding.rf.index, 4] <- same.item.in.a.row.counter
      same.item.in.a.row.counter <- NULL
      if (!is.null(holdout.frequencies)) {
        RF.matrix[corresponding.rf.index, 5] <- x.star.total
        x.star.total = NULL
      }
    }
  }
  return(RF.matrix)
}

