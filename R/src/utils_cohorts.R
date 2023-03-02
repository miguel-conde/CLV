
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
    mutate(date = as_date(date)) %>% 
    mutate(sales = ifelse(is.na(sales), 0, sales))
  
  return(out)
}

build_rf_table <- function(probe) {
  
  aux <- probe %>% spread(date, sales) %>% select(-id) 
  
  out <- tibble(f = aux %>% rowSums(na.rm = TRUE)-1,
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

my_bgbb.PlotFreqVsConditionalExpectedFrequency <- 
  function (params, n.star, rf.matrix, x.star, trunc = NULL, xlab = "Calibration period transactions", 
            ylab = "Holdout period transactions", xticklab = NULL, title = "Conditional Expectation") 
  {
  if (length(x.star) != nrow(rf.matrix)) 
    stop("x.star must have the same number of entries as rows in rf.matrix")
  if (!(length(n.star) == 1 || length(n.star) == nrow(rf.matrix))) 
    stop("n.star must be a single value or have as many entries as rows in rf.matrix")
  dc.check.model.params(printnames = c("r", "alpha", "s", "beta"), 
                        params = params, func = "bgbb.PlotFreqVsConditionalExpectedFrequency")
  if (any(x.star < 0) || !is.numeric(x.star)) 
    stop("x.star must be numeric and may not contain negative numbers.")
  if (any(n.star < 0) || !is.numeric(n.star)) 
    stop("n.star must be numeric and may not contain negative numbers.")
  n.star <- rep(n.star, length.out = nrow(rf.matrix))
  tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a frequency column labelled \"x\""))
  tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a recency column labelled \"t.x\""))
  tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
  tryCatch(custs <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a column for the number of customers that have each combination of \"x\", \"t.x\", and \"n.cal\", labelled \"custs\""))
  if (is.null(trunc)) 
    trunc <- max(n.cal)
  if (trunc > max(n.cal)) {
    warning("The truncation number provided in bgbb.PlotFreqVsConditionalExpectedFrequency was greater than the maximum number of possible transactions. It has been reduced to ", 
            max(n.cal))
    trunc = max(n.cal)
  }
  actual.freq <- rep(0, max(n.cal) + 1)
  expected.freq <- rep(0, max(n.cal) + 1)
  bin.size <- rep(0, max(n.cal) + 1)
  for (ii in 0:max(n.cal)) {
    bin.size[ii + 1] <- sum(custs[x == ii])
    actual.freq[ii + 1] <- sum(x.star[x == ii])
    expected.freq[ii + 1] <- sum(bgbb.ConditionalExpectedTransactions(params, 
                                                                      n.cal[x == ii], n.star[x == ii], ii, t.x[x == ii]) * 
                                   custs[x == ii])
  }
  comparison <- rbind(actual.freq/bin.size, expected.freq/bin.size, 
                      bin.size)
  colnames(comparison) <- paste("freq.", 0:max(n.cal), sep = "")
  if (is.null(xticklab) == FALSE) {
    x.labels <- xticklab
  }
  else {
    if ((trunc + 1) < ncol(comparison)) {
      x.labels <- 0:(trunc)
    }
    else {
      x.labels <- 0:(ncol(comparison) - 1)
    }
  }
  actual.freq <- comparison[1, 1:(trunc + 1)]
  expected.freq <- comparison[2, 1:(trunc + 1)]
  custs.in.plot <- sum(comparison[3, 1:(trunc + 1)])
  if (custs.in.plot < 0.9 * sum(custs)) {
    warning("Less than 90% of customers are represented in your plot (", 
            custs.in.plot, " of ", sum(custs), " are plotted).")
  }
  # Faltaba na.rm = TRUE
  ylim <- c(0, ceiling(max(c(actual.freq, expected.freq), na.rm = TRUE) * 
                         1.1))
  plot(actual.freq, type = "l", xaxt = "n", col = 1, ylim = ylim, 
       xlab = xlab, ylab = ylab, main = title)
  lines(expected.freq, lty = 2, col = 2)
  axis(1, at = 1:(trunc + 1), labels = x.labels)
  legend("topleft", legend = c("Actual", "Model"), col = 1:2, 
         lty = 1:2, lwd = 1)
  return(comparison)
}

my_bgbb.PlotRecVsConditionalExpectedFrequency <- 
  function (params, n.star, rf.matrix, x.star, trunc = NULL, xlab = "Calibration period recency", 
            ylab = "Holdout period transactions", xticklab = NULL, title = "Conditional Expected Transactions by Recency") 
  {
  if (length(x.star) != nrow(rf.matrix)) 
    stop("x.star must have the same number of entries as rows in rf.matrix")
  if (!(length(n.star) == 1 || length(n.star) == nrow(rf.matrix))) 
    stop("n.star must be a single value or have as many entries as rows in rf.matrix")
  dc.check.model.params(printnames = c("r", "alpha", "s", "beta"), 
                        params = params, func = "bgbb.PlotRecVsConditionalExpectedFrequency")
  if (any(x.star < 0) || !is.numeric(x.star)) 
    stop("x.star must be numeric and may not contain negative numbers.")
  if (any(n.star < 0) || !is.numeric(n.star)) 
    stop("n.star must be numeric and may not contain negative numbers.")
  n.star <- rep(n.star, length.out = nrow(rf.matrix))
  tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a frequency column labelled \"x\""))
  tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a recency column labelled \"t.x\""))
  tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
  tryCatch(custs <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a column for the number of customers that have each combination of \"x\", \"t.x\", and \"n.cal\", labelled \"custs\""))
  if (is.null(trunc)) 
    trunc <- max(n.cal)
  if (trunc > max(n.cal)) {
    warning("The truncation number provided in bgbb.PlotRecVsConditionalExpectedFrequency was greater than the maximum number of possible transactions. It has been reduced to ", 
            max(n.cal))
    trunc = max(n.cal)
  }
  actual.freq <- rep(0, max(n.cal) + 1)
  expected.freq <- rep(0, max(n.cal) + 1)
  bin.size <- rep(0, max(n.cal) + 1)
  for (ii in 0:max(n.cal)) {
    bin.size[ii + 1] <- sum(custs[t.x == ii])
    actual.freq[ii + 1] <- sum(x.star[t.x == ii])
    expected.freq[ii + 1] <- sum(bgbb.ConditionalExpectedTransactions(params, 
                                                                      n.cal[t.x == ii], n.star[t.x == ii], x[t.x == ii], 
                                                                      ii) * custs[t.x == ii])
  }
  comparison <- rbind(actual.freq/bin.size, expected.freq/bin.size, 
                      bin.size)
  colnames(comparison) <- paste("rec.", 0:max(n.cal), sep = "")
  custs.in.plot <- sum(comparison[3, 1:(trunc + 1)])
  if (custs.in.plot < 0.9 * sum(custs)) {
    warning("Less than 90% of customers are represented in your plot (", 
            custs.in.plot, " of ", sum(custs), " are plotted).")
  }
  if (is.null(xticklab) == FALSE) {
    x.labels <- xticklab
  }
  else {
    if ((trunc + 1) < ncol(comparison)) {
      x.labels <- 0:(trunc)
    }
    else {
      x.labels <- 0:(ncol(comparison) - 1)
    }
  }
  actual.freq <- comparison[1, 1:(trunc + 1)]
  expected.freq <- comparison[2, 1:(trunc + 1)]
  # Faltaba na.rm = TRUE
  ylim <- c(0, ceiling(max(c(actual.freq, expected.freq), na.rm = TRUE) * 
                         1.1))
  plot(actual.freq, type = "l", xaxt = "n", col = 1, ylim = ylim, 
       xlab = xlab, ylab = ylab, main = title)
  lines(expected.freq, lty = 2, col = 2)
  axis(1, at = 1:(trunc + 1), labels = x.labels)
  legend("topleft", legend = c("Actual", "Model"), col = 1:2, 
         lty = 1:2, lwd = 1)
  return(comparison)
}

####
# my_dc.ElogToCbsCbt <- function (elog, per = "week", T.cal = max(elog$date), T.tot = max(elog$date), 
#                                 merge.same.date = TRUE, cohort.birth.per = T.cal, dissipate.factor = 1, 
#                                 statistic = "freq") 
# {
#   dc.WriteLine("Started making CBS and CBT from the ELOG...")
#   elog <- dc.FilterCustByBirth(elog, cohort.birth.per)
#   if (nrow(elog) == 0) 
#     stop("error caused by customer birth filtering")
#   elog <- elog[elog$date <= T.tot, ]
#   if (nrow(elog) == 0) 
#     stop("error caused by holdout period end date")
#   elog <- dc.DissipateElog(elog, dissipate.factor)
#   if (nrow(elog) == 0) 
#     stop("error caused by event long dissipation")
#   if (merge.same.date) {
#     elog <- dc.MergeTransactionsOnSameDate(elog)
#     if (nrow(elog) == 0) 
#       stop("error caused by event log merging")
#   }
#   calibration.elog <- elog[elog$date <= T.cal, ]
#   holdout.elog <- elog[elog$date > T.cal, ]
#   split.elog.list <- dc.SplitUpElogForRepeatTrans(calibration.elog)
#   repeat.transactions.elog <- split.elog.list$repeat.trans.elog
#   cust.data <- split.elog.list$cust.data
#   dc.WriteLine("Started Building CBS and CBT for calibration period...")
#   cbt.cal <- dc.BuildCBTFromElog(calibration.elog, statistic)
#   cbt.cal.rep.trans <- dc.BuildCBTFromElog(repeat.transactions.elog, 
#                                            statistic)
#   cbt.cal <- dc.MergeCustomers(cbt.cal, cbt.cal.rep.trans)
#   dates <- data.frame(cust.data$birth.per, cust.data$last.date, 
#                       T.cal)
#   cbs.cal <- my_dc.BuildCBSFromCBTAndDates(cbt.cal, dates, per, 
#                                            cbt.is.during.cal.period = TRUE)
#   dc.WriteLine("Finished building CBS and CBT for calibration period.")
#   cbt.holdout <- NULL
#   cbs.holdout <- NULL
#   if (nrow(holdout.elog) > 0) {
#     dc.WriteLine("Started building CBS and CBT for holdout period...")
#     cbt.holdout <- dc.BuildCBTFromElog(holdout.elog, statistic)
#     dates <- c((T.cal + 1), T.tot)
#     cbs.holdout <- my_dc.BuildCBSFromCBTAndDates(cbt.holdout, 
#                                                  dates, per, cbt.is.during.cal.period = FALSE)
#     cbt.holdout <- dc.MergeCustomers(cbt.cal, cbt.holdout)
#     cbs.holdout <- dc.MergeCustomers(cbs.cal, cbs.holdout)
#     dc.WriteLine("Finished building CBS and CBT for holdout.")
#     dc.WriteLine("...Finished Making All CBS and CBT")
#     return(list(cal = list(cbs = cbs.cal, cbt = cbt.cal), 
#                 holdout = list(cbt = cbt.holdout, cbs = cbs.holdout), 
#                 cust.data = cust.data))
#   }
#   dc.WriteLine("...Finished Making All CBS and CBT")
#   return(list(cal = list(cbs = cbs.cal, cbt = cbt.cal), holdout = list(cbt = cbt.holdout, 
#                                                                        cbs = cbs.holdout), cust.data = cust.data))
# }
# 
# my_dc.BuildCBSFromCBTAndDates <- function (cbt, dates, per, cbt.is.during.cal.period = TRUE) 
# {
#   if (cbt.is.during.cal.period == TRUE) {
#     dc.WriteLine("Started making calibration period CBS...")
#     custs.first.dates <- dates[, 1]
#     custs.last.dates <- dates[, 2]
#     T.cal <- dates[, 3]
#     if (length(custs.first.dates) != length(custs.last.dates)) {
#       stop("Invalid dates (different lengths) in BuildCBSFromFreqCBTAndDates")
#     }
#     f <- rowSums(cbt)
#     r <- as.numeric(difftime(custs.last.dates, custs.first.dates, 
#                              units = "days"))
#     T <- as.numeric(difftime(T.cal, custs.first.dates, units = "days"))
#     x <- switch(per, day = 1, week = 7, month = 365/12, quarter = 365/4, 
#                 year = 365)
#     r = r/x
#     T = T/x
#     cbs = cbind(f, r, T)
#     rownames(cbs) <- rownames(cbt)
#     colnames(cbs) <- c("x", "t.x", "T.cal")
#   }
#   else {
#     dc.WriteLine("Started making holdout period CBS...")
#     date.begin.holdout.period <- dates[1]
#     date.end.holdout.period <- dates[2]
#     f <- rowSums(cbt)
#     T <- as.numeric(difftime(date.end.holdout.period, date.begin.holdout.period, 
#                              units = "days")) + 1
#     x <- switch(per, day = 1, week = 7, month = 365/12, quarter = 365/4, 
#                 year = 365)
#     T = T/x
#     cbs = cbind(f, T)
#     rownames(cbs) <- rownames(cbt)
#     colnames(cbs) <- c("x.star", "T.star")
#   }
#   dc.WriteLine("Finished building CBS.")
#   return(cbs)
# }
