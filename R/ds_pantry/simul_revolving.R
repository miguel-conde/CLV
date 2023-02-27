library(tidyverse)
library(lubridate)

source(here::here("R", "src", "utils_cohorts.R"), encoding = "utf8")



curve(dbeta(x, .1, .9), 0, 1)
# - lambda - exponencial(lambda = 1/5)
rpois(10, rexp(10, 1/5))




simul_feb <- simul_N("2020-02-01", "2022-12-31", N = 1000)
rf_feb <- simul_feb %>% build_rf_table() 

# Event log
elog_feb <- simul_feb %>% 
  # rename(cust = id) %>% 
  mutate(sales = ifelse(is.na(sales), 0, sales)) %>% 
  filter(sales != 0)

rf_matrix_feb <- elog_feb %>% 
  build_rf_table()  %>% 
  rename(x = f, t.x = r, custs = n, n.cal = T) %>% 
  as.matrix()

# STEP BY STEP ------------------------------------------------------------

# Event log
elog <- elog_feb %>% rename(cust = id)
elog

# Merge all transactions that occurred on the same day
elog <- dc.MergeTransactionsOnSameDate(elog)

# Calibration & holdout
end.of.cal.period <- as.Date("2021-12-31")
elog.cal <- elog[which(elog$date <= end.of.cal.period), ]

# the model is generally concerned with repeat transactions—that
# is, the first transaction is ignored.
# The one problem with simply getting rid of customers’ first transactions is 
# the following: We have to keep track of a “time zero” as a point of reference 
# for recency and total time observed. For this
# reason, we use dc.SplitUpElogForRepeatTrans, which returns a filtered event
# log ($repeat.trans.elog) as well as saving important information about each
# customer ($cust.data).



# Intento 1 ---------------------------------------------------------------



library(BTYD)

## Data Preparation

rf_matrix <- rf_feb %>% 
  rename(x = f, t.x = r, n.cal = `T`, custs = n)

## Parameter Estimation

params   <- bgbb.EstimateParameters(rf.matrix)
LL       <- bgbb.rf.matrix.LL(params, rf.matrix)
p.matrix <- c(params, LL)

for (i in 1:2){
  params       <- bgbb.EstimateParameters(rf.matrix, params);
  LL           <- bgbb.rf.matrix.LL(params, rf.matrix);
  p.matrix.row <- c(params, LL);
  p.matrix     <- rbind(p.matrix, p.matrix.row);
}
colnames(p.matrix) <- c("alpha", "beta", "gamma", "delta", "LL")
rownames(p.matrix) <- 1:3
p.matrix

bgbb.PlotTransactionRateHeterogeneity(params)
bgbb.PlotDropoutRateHeterogeneity(params)

## Individual Level Estimations

# number of repeat transactions we expect a newly-acquired customer to make in 
# a period of 36 months
bgbb.Expectation(params, n = 36)

# A, who made 0 transactions in the calibration period; 
# customer A
n.cal  = 7  # number of transaction opportunities in the calibration period.
n.star = 10 # number of transaction opportunities in the holdout period,
x      = 0  # number of repeat transactions made by the customer in the calibration period.
t.x    = 0  # ecency - the transaction opportunity in which the customer made their last transaction
bgbb.ConditionalExpectedTransactions(params, n.cal,
                                     n.star, x, t.x)


# and B, who made 4 transactions in the calibration period, with the last 
# transaction occuring in the 5th month.
# customer B
x   = 23
t.x = 30
bgbb.ConditionalExpectedTransactions(params, n.cal,
                                     n.star, x, t.x)
## Plotting/ Goodness-of-fit

bgbb.PlotFrequencyInCalibration(params, rf.matrix)



# Intento 2 ---------------------------------------------------------------

# Event log
elog_feb <- simul_feb %>% 
  rename(cust = id) %>% 
  mutate(sales = coalesce(is.na(sales), 0)) %>% 
  filter(sales != 0)

cal.rf.matrix <- rf_feb %>% 
  rename(x = f, t.x = r, custs = n, n.cal = T) %>% 
  as.matrix()

max(elog_feb$date);
min(elog_feb$date);

T.cal <- as.Date("2021-12-31")
simData <- dc.ElogToCbsCbt(elog_feb, per="month", T.cal)

# CALIBRATION RF MATRIX
cal.cbs <- simData$cal$cbs
freq    <- cal.cbs[,"x"]
rec     <- cal.cbs[,"t.x"]

trans.opp <- 24 # transaction opportunities
cal.rf.matrix <- dc.MakeRFmatrixCal(freq, rec, trans.opp)
cal.rf.matrix[1:5,]

# HOLDOUT RF MATRIX
holdout.rf.matrix <- dc.MakeRFmatrixHoldout(simData$holdout$cbt)

##

params   <- bgbb.EstimateParameters(cal.rf.matrix)
LL       <- bgbb.rf.matrix.LL(params, cal.rf.matrix)
p.matrix <- c(params, LL)

for (i in 1:2){
  params       <- bgbb.EstimateParameters(cal.rf.matrix, params)
  LL           <- bgbb.rf.matrix.LL(params, cal.rf.matrix)
  p.matrix.row <- c(params, LL)
  p.matrix     <- rbind(p.matrix, p.matrix.row)
}

colnames(p.matrix) <- c("alpha", "beta", "gamma", "delta", "LL")
rownames(p.matrix) <- 1:3
p.matrix

bgbb.PlotTransactionRateHeterogeneity(params)
bgbb.PlotDropoutRateHeterogeneity(params)

## Individual Level Estimations

# number of repeat transactions we expect a newly-acquired customer to make in 
# a period of 36 months
bgbb.Expectation(params, n = 36)

# A, who made 0 transactions in the calibration period; 
# customer A
n.cal  = 24 # number of transaction opportunities in the calibration period.
n.star = 12 # number of transaction opportunities in the holdout period,
x      = 0  # number of repeat transactions made by the customer in the calibration period.
t.x    = 0  # recency - the transaction opportunity in which the customer made their last transaction
bgbb.ConditionalExpectedTransactions(params, n.cal,
                                     n.star, x, t.x)


# and B, who made 4 transactions in the calibration period, with the last 
# transaction occuring in the 5th month.
# customer B
x   = 23
t.x = 30
bgbb.ConditionalExpectedTransactions(params, n.cal,
                                     n.star, x, t.x)
## Plotting/ Goodness-of-fit

bgbb.PlotFrequencyInCalibration(params, cal.rf.matrix)

cal.rf.matrix %>% as_tibble() %>% group_by(x) %>% summarise(custs = sum(custs)) %>% as.data.frame()

####
# holdout.cbs <- simData$holdout$cbs
# x.star <- holdout.cbs[,"x.star"]

# holdout.cbs <- simData$holdout$cbs
# freq    <- holdout.cbs[,"x"]
# rec     <- holdout.cbs[,"t.x"]
# 
# trans.opp <- 12 # transaction opportunities
# holdout.cbs.rf.matrix <- dc.MakeRFmatrixCal(freq, rec, trans.opp)

# Actual vs. conditional expected transactions in the holdout period,
# binned by calibration period frequency
n.star <- 12 # length of the holdout period
x.star <- holdout.rf.matrix[, "x.star"]
comp <- bgbb.PlotFreqVsConditionalExpectedFrequency(params, n.star,
                                                    cal.rf.matrix, x.star)

rownames(comp) <- c("act", "exp", "bin")
comp

# Actual vs. conditional expected transactions in the holdout period,
# binned by calibration period recency.
comp <- bgbb.PlotRecVsConditionalExpectedFrequency(params, n.star,
                                                   rf.matrix, x.star)
rownames(comp) <- c("act", "exp", "bin")
comp

# Actual vs. expected incremental purchasing behaviour.
inc.track.data <- donationsSummary$annual.trans
n.cal <- 6
xtickmarks <- 1996:2006
inc.tracking <- bgbb.PlotTrackingInc(params, rf.matrix,
                                     inc.track.data,
                                     xticklab = xtickmarks)
rownames(inc.tracking) <- c("act", "exp")
inc.tracking

# Actual vs. expected cumulative purchasing behaviour.
cum.track.data <- cumsum(inc.track.data)
cum.tracking <- bgbb.PlotTrackingCum(params, rf.matrix, cum.track.data,
                                     xticklab = xtickmarks)

rownames(cum.tracking) <- c("act", "exp")
cum.tracking

