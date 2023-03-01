library(tidyverse)
library(lubridate)
library(BTYD)

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

# Data Preparation --------------------------------------------------------


# Event log
elog <- elog_feb %>% rename(cust = id)
elog

# Merge all transactions that occurred on the same day
elog <- dc.MergeTransactionsOnSameDate(elog)
elog %>% head()

# Calibration & holdout
end.of.cal.period <- as.Date("2021-12-31")
elog.cal <- elog[which(elog$date <= end.of.cal.period), ]
elog.cal %>% head()

# the model is generally concerned with repeat transactions—that
# is, the first transaction is ignored.
# The one problem with simply getting rid of customers’ first transactions is 
# the following: We have to keep track of a “time zero” as a point of reference 
# for recency and total time observed. For this
# reason, we use dc.SplitUpElogForRepeatTrans, which returns a filtered event
# log ($repeat.trans.elog) as well as saving important information about each
# customer ($cust.data).
split.data <- dc.SplitUpElogForRepeatTrans(elog.cal)

split.data$repeat.trans.elog %>% head() # filtered event log, the first transaction is ignored
split.data$cust.data %>% head() # mportant information about each customer

clean.elog <- split.data$repeat.trans.elog

# The next step is to create a Customer-By-Time matrix (CBT)
# If you have already used dc.MergeTransactionsOnSameDate, this will simply be a
# reach customer-by-time matrix.
freq.cbt <- dc.CreateFreqCBT(clean.elog)
freq.cbt[1:3,1:7]

dc.CreateReachCBT(clean.elog)[1:3,1:7] # Contains a 1 if the customer made any
                                       # transactions on that day, and 0 otherwise.
dc.CreateSpendCBT(clean.elog)[1:3,1:7] # Contains the amount spent by that customer
                                       # on that day.

# We have a small problem now—since we have deleted all the first transactions,
# the frequency customer-by-time matrix does not have any of the customers who
# made zero repeat transactions. These customers are still important; in fact, in
# most datasets, more customers make zero repeat transactions than any other
# number. Solving the problem is reasonably simple: we create a customer-by-time
# matrix using all transactions, and then merge the filtered CBT with this total
# CBT (using data from the filtered CBT and customer IDs from the total CBT).
tot.cbt <- dc.CreateFreqCBT(elog)
tot.cbt[1:3,1:7]
cal.cbt <- dc.MergeCustomers(tot.cbt, freq.cbt) # Takes two CBT or CBS matrices 
                                                # and ensures that the second one 
                                                # has the same row names as the first.
cal.cbt[1:3,1:7]
plot(as_date(colnames(cal.cbt)), cal.cbt %>% colSums() )

# En nuestro caso ningún cliente tiene solo 1  transaccion
freq.cbt %>% rowSums() %>% is.zero(.) %>% sum()
tot.cbt %>% rowSums() %>% sum() %>% is.zero(.) %>% sum()
cal.cbt %>% rowSums() %>% sum() %>% is.zero(.) %>% sum()

elog %>% group_by(cust) %>% summarise(n = sum(sales)) %>% arrange(n)

# Customer-By-Sufficient-statistic matrix (CBS)
#  The function we are going to use is dc.BuildCBSFromCBTAndDates, which 
# requires:
#   - a customer-by-time matrix,
#   - starting and ending dates for each customer, 
#   - and the time of the end of the calibration period. 
#   - It also requires a time period to use—in this case, we are choosing to use 
#     weeks. 
#
# This function could also be used for the holdout period by setting 
# cbt.is.during.cal.period to FALSE. There are slight differences when this
# function is used for the holdout period—it requires different input dates (simply
# the start and end of the holdout period) and does not return a recency (which
# has little value in the holdout period).
birth.periods <- split.data$cust.data$birth.per
last.dates    <- split.data$cust.data$last.date
cal.cbs.dates <- data.frame(birth.periods, last.dates, end.of.cal.period)

# The customer-by-sufficient statistic matrix will contain:
#   - the sum of the statistic included in the customer-by-time matrix,
#   - the customer's last transaction date, 
#   - and the total time period for which the customer was observed. 
#
# Customer-by-sufficient-statistic matrix (CBS), with three columns: 
#   - frequency("x"), 
#   - recency("t.x") 
#   - and total time observed("T.cal"). 
# Frequency is total transactions, not repeat transactions.
cal.cbs       <- dc.BuildCBSFromCBTAndDates(cal.cbt, cal.cbs.dates, per="month")

cal.cbs %>% head()

# Create a recencyfrequency matrix from our customer-by-sufficient-statistic matrix
freq<- cal.cbs[,"x"]
rec <- cal.cbs[,"t.x"]
trans.opp <- 22 # transaction opportunities, cohorte feb => 23 - 1 (1ª transacción)
cal.rf.matrix <- my_dc.MakeRFmatrixCal(freq, rec, trans.opp)

# Parameter Estimation ----------------------------------------------------

params <- bgbb.EstimateParameters(cal.rf.matrix)

LL <- bgbb.rf.matrix.LL(params, cal.rf.matrix)
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


# Individual level estimation ---------------------------------------------

# number of transactions we expect a newly acquired customer
# to make in a given time period
bgbb.Expectation(params, n=10)

# say something about our existing customers, not
# just about a hypothetical customer to be acquired in the future

# A, who made 0 transactions in the calibration period; 
# customer A
n.cal  = 22 # number of transaction opportunities in the calibration period.
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

# We can also obtain expected characteristics for a specific customer, conditional 
# on their purchasing behavior during the calibration period.
cal.cbs["516",]

x     <- cal.cbs["516", "x"]
t.x   <- cal.cbs["516", "t.x"]
n.cal <- cal.cbs["516", "T.cal"]
bgbb.ConditionalExpectedTransactions(params, n.cal,
                                     n.star, x, t.x)

bgbb.PAlive(params, x, t.x, n.cal)


# Plotting/ Goodness-of-fit -----------------------------------------------

bgbb.PlotFrequencyInCalibration(params = params,
                                rf.matrix = cal.rf.matrix, 
                                censor = 15)

## Holdout
# dc.ElogToCbsCbt
# produces both a calibration period customer-by-sufficient-statistic matrix and a
# holdout period customer-by-sufficient-statistic matrix, which could be combined
# in order to find the number of transactions each customer made in the holdout
# period.
#  Note that I subtract the number of repeat
# transactions in the calibration period from the total number of transactions. We
# remove the initial transactions first as we are not concerned with them.

#####

# Auxiliary Event log - nos quedamos solo con las trasacciones de repetición
# de TODO el event lo original
aux_elog <- dc.SplitUpElogForRepeatTrans(elog)$repeat.trans.elog

# En x.start vamos a añadir a cal.cbs las frecuencias de holdout
x.star <- rep(0, nrow(cal.cbs))
cal.cbs <- cbind(cal.cbs, x.star)
elog.custs <- aux_elog$cust

# Recorremos por fila cal.cbs
for (i in 1:nrow(cal.cbs)){
  # Para el cliente i de la CBS calibrado...
  current.cust         <- rownames(cal.cbs)[i]
  # ...miro el total de transacciones en calibrado + holdout,
  tot.cust.trans       <- length(which(elog.custs == current.cust))
  # miro también las transacciones solo de calibrado...
  cal.trans            <- cal.cbs[i, "x"]
  # ... en la columna x.start de la CBS de calibrado meto la diferencia,
  # que son las transacciones de holdout
  cal.cbs[i, "x.star"] <- tot.cust.trans - cal.trans
}

round(cal.cbs, digits = 3)

#####
cal_holdout_data <- dc.ElogToCbsCbt(elog, per="month", end.of.cal.period)

cal.rf.matrix <- 
  my_dc.MakeRFmatrixCal(frequencies = cal_holdout_data$cal$cbs[,"x"], 
                        periods.of.final.purchases = cal_holdout_data$cal$cbs[,"t.x"], 
                        num.of.purchase.periods = trans.opp)

holdout.rf.matrix <- 
  my_dc.MakeRFmatrixCal(frequencies = cal_holdout_data$cal$cbs[,"x"], 
                        periods.of.final.purchases = cal_holdout_data$cal$cbs[,"t.x"], 
                        num.of.purchase.periods = trans.opp, 
                        holdout.frequencies = cal_holdout_data$holdout$cbs[,"x.star"])

# Now we can see how well our model does in the holdout period.
comp <- my_bgbb.PlotFreqVsConditionalExpectedFrequency(params = params, 
                                                       n.star = 12, 
                                                       rf.matrix = cal.rf.matrix,
                                                       x.star = holdout.rf.matrix[, "x.star"])

rownames(comp) <- c("act", "exp", "bin")
round(comp, digits = 3)

##
comp <- my_bgbb.PlotRecVsConditionalExpectedFrequency(params = params, 
                                                      n.star = 12, 
                                                      rf.matrix = cal.rf.matrix,
                                                      x.star = holdout.rf.matrix[, "x.star"])
rownames(comp) <- c("act", "exp", "bin")
comp

##
inc.track.data <- elog %>% group_by(year(date), month(date)) %>% summarise(n = sum(sales)) %>% pull(n)
# xtickmarks <- 1996:2006
xtickmarks <- 1:35
inc.tracking <- bgbb.PlotTrackingInc(params = params, 
                                     rf.matrix = cal.rf.matrix,
                                     actual.inc.repeat.transactions = inc.track.data,
                                     xticklab = xtickmarks)

rownames(inc.tracking) <- c("act", "exp")
inc.tracking

##
cum.track.data <- cumsum(inc.track.data)

cum.tracking <- bgbb.PlotTrackingCum(params                         = params, 
                                     rf.matrix                      = cal.rf.matrix, 
                                     actual.cum.repeat.transactions = cum.track.data,
                                     xticklab                       = xtickmarks)

rownames(cum.tracking) <- c("act", "exp")
cum.tracking


# AL GRANO ----------------------------------------------------------------

# Data preparation --------------------------------------------------------

# Event log
elog <- elog_feb %>% rename(cust = id)
elog

end.of.cal.period <- as.Date("2021-12-31")

cal_holdout_data <- dc.ElogToCbsCbt(elog, per="month", end.of.cal.period)
trans.opp <- 22 # transaction opportunities, cohorte feb => 23 - 1 (1ª transacción)

rf.matrix <- 
  my_dc.MakeRFmatrixCal(frequencies                = cal_holdout_data$cal$cbs[,"x"], 
                        periods.of.final.purchases = cal_holdout_data$cal$cbs[,"t.x"], 
                        num.of.purchase.periods    = trans.opp, 
                        holdout.frequencies        = cal_holdout_data$holdout$cbs[,"x.star"])

# Parameter Estimation ----------------------------------------------------

params <- bgbb.EstimateParameters(cal.rf.matrix)

LL <- bgbb.rf.matrix.LL(params, cal.rf.matrix)
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


# Individual level estimation ---------------------------------------------

# number of transactions we expect a newly acquired customer
# to make in a given time period
bgbb.Expectation(params, n=10)

# say something about our existing customers, not
# just about a hypothetical customer to be acquired in the future

# A, who made 0 transactions in the calibration period; 
# customer A
n.cal  = 22 # number of transaction opportunities in the calibration period.
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

# We can also obtain expected characteristics for a specific customer, conditional 
# on their purchasing behavior during the calibration period.
cal.cbs["516",]

x     <- cal.cbs["516", "x"]
t.x   <- cal.cbs["516", "t.x"]
n.cal <- cal.cbs["516", "T.cal"]
bgbb.ConditionalExpectedTransactions(params, n.cal,
                                     n.star, x, t.x)

bgbb.PAlive(params, x, t.x, n.cal)

# Plotting/ Goodness-of-fit -----------------------------------------------

bgbb.PlotFrequencyInCalibration(params = params,
                                rf.matrix = cal.rf.matrix, 
                                censor = 15)

# Now we can see how well our model does in the holdout period.
comp <- my_bgbb.PlotFreqVsConditionalExpectedFrequency(params = params, 
                                                       n.star = 12, 
                                                       rf.matrix = rf.matrix,
                                                       x.star = rf.matrix[, "x.star"])

rownames(comp) <- c("act", "exp", "bin")
round(comp, digits = 3)

##
comp <- my_bgbb.PlotRecVsConditionalExpectedFrequency(params = params, 
                                                      n.star = 12, 
                                                      rf.matrix = rf.matrix,
                                                      x.star = rf.matrix[, "x.star"])
rownames(comp) <- c("act", "exp", "bin")
comp

##
inc.track.data <- elog %>% group_by(year(date), month(date)) %>% summarise(n = sum(sales)) %>% pull(n)
xtickmarks <- 1:35
inc.tracking <- bgbb.PlotTrackingInc(params = params, 
                                     rf.matrix = rf.matrix,
                                     actual.inc.repeat.transactions = inc.track.data,
                                     xticklab = xtickmarks)

rownames(inc.tracking) <- c("act", "exp")
inc.tracking

##
cum.track.data <- cumsum(inc.track.data)

cum.tracking <- bgbb.PlotTrackingCum(params                         = params, 
                                     rf.matrix                      = rf.matrix, 
                                     actual.cum.repeat.transactions = cum.track.data,
                                     xticklab                       = xtickmarks)

rownames(cum.tracking) <- c("act", "exp")
cum.tracking


