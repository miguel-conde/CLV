library(tidyverse)
library(BTYD)


# BG/BB -------------------------------------------------------------------

# A recency-frequency matrix consists of a
# row for every possible calibration period recency/frequency combination, and
# contains the total number of customers which had that particular combination
# of recency and frequency.

simElog <- system.file("data/discreteSimElog.csv",
                       package = "BTYD")
elog <- dc.ReadLines(simElog, cust.idx = 1, date.idx = 2)
elog[1:3,]
# cust date
# 1 1 1970-01-01
# 2 1 1975-01-01
# 3 1 1977-01-01
elog$date <- as.Date(elog$date, "%Y-%m-%d")
max(elog$date);
# [1] "1983-01-01"
min(elog$date);
# [1] "1970-01-01"
# let’s make the calibration period end somewhere in-between
T.cal <- as.Date("1977-01-01")
simData <- dc.ElogToCbsCbt(elog, per="year", T.cal)
cal.cbs <- simData$cal$cbs
freq<- cal.cbs[,"x"]
rec <- cal.cbs[,"t.x"]
trans.opp <- 7 # transaction opportunities
cal.rf.matrix <- dc.MakeRFmatrixCal(freq, rec, trans.opp)
cal.rf.matrix[1:5,]
# x t.x n.cal custs
# [1,] 0 0 7 2900
# [2,] 1 1 7 933
# [3,] 1 2 7 218
# [4,] 2 2 7 489
# [5,] 1 3 7 95

data(donationsSummary)

rf.matrix <- donationsSummary$rf.matrix
params <- bgbb.EstimateParameters(rf.matrix);
LL <- bgbb.rf.matrix.LL(params, rf.matrix);
p.matrix <- c(params, LL);
for (i in 1:2){
  params <- bgbb.EstimateParameters(rf.matrix, params);
  LL <- bgbb.rf.matrix.LL(params, rf.matrix);
  p.matrix.row <- c(params, LL);
  p.matrix <- rbind(p.matrix, p.matrix.row);
}
colnames(p.matrix) <- c("alpha", "beta", "gamma", "delta", "LL");
rownames(p.matrix) <- 1:3;
p.matrix;

bgbb.Expectation(params, n=10);

# customer A
n.cal = 6
n.star = 10
x = 0
t.x = 0
bgbb.ConditionalExpectedTransactions(params, n.cal,
                                     n.star, x, t.x)
# [1] 0.1302169
# customer B
x = 4
t.x = 5
bgbb.ConditionalExpectedTransactions(params, n.cal,
                                     n.star, x, t.x)
# [1] 3.627858

bgbb.PlotFrequencyInCalibration(params, rf.matrix)

holdout.cbs <- simData$holdout$cbs
x.star <- holdout.cbs[,"x.star"]

# Actual vs. conditional expected transactions in the holdout period,
# binned by calibration period frequency
n.star <- 5 # length of the holdout period
x.star <- donationsSummary$x.star
comp <- bgbb.PlotFreqVsConditionalExpectedFrequency(params, n.star,
                                                    rf.matrix, x.star)

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


# Pareto/NBD --------------------------------------------------------------

cdnowElog <- system.file("data/cdnowElog.csv", package = "BTYD")
elog <- dc.ReadLines(cdnowElog, cust.idx = 2,
                     date.idx = 3, sales.idx = 5)
elog[1:3,]


elog$date <- as.Date(elog$date, "%Y%m%d")
elog[1:3,]

elog <- dc.MergeTransactionsOnSameDate(elog)

end.of.cal.period <- as.Date("1997-09-30")
elog.cal <- elog[which(elog$date <= end.of.cal.period), ]

split.data <- dc.SplitUpElogForRepeatTrans(elog.cal)
clean.elog <- split.data$repeat.trans.elog

freq.cbt <- dc.CreateFreqCBT(clean.elog)
freq.cbt[1:3,1:5]

tot.cbt <- dc.CreateFreqCBT(elog)
cal.cbt <- dc.MergeCustomers(tot.cbt, freq.cbt)

birth.periods <- split.data$cust.data$birth.per
last.dates    <- split.data$cust.data$last.date
cal.cbs.dates <- data.frame(birth.periods, last.dates, end.of.cal.period)
cal.cbs       <- dc.BuildCBSFromCBTAndDates(cal.cbt, cal.cbs.dates, per="week")

params <- pnbd.EstimateParameters(cal.cbs = cal.cbs)
round(params, digits = 3)
# [1] 0.553 10.580 0.606 11.656
LL <- pnbd.cbs.LL(params = params,
                  cal.cbs = cal.cbs)
LL


p.matrix <- c(params, LL)
for (i in 1:2){
  params <- pnbd.EstimateParameters(cal.cbs = cal.cbs,
                                    par.start = params)
  LL <- pnbd.cbs.LL(params = params,
                    cal.cbs = cal.cbs)
  p.matrix.row <- c(params, LL)
  p.matrix <- rbind(p.matrix, p.matrix.row)
}
colnames(p.matrix) <- c("r", "alpha", "s", "beta", "LL")
rownames(p.matrix) <- 1:3
round(p.matrix, digits = 3)

pnbd.PlotTransactionRateHeterogeneity(params)
pnbd.PlotDropoutRateHeterogeneity(params)

pnbd.Expectation(params = params, t = 52)

cal.cbs["1516",]
# x t.x T.cal
# 26.00000 30.85714 31.00000
x <- cal.cbs["1516", "x"]
t.x <- cal.cbs["1516", "t.x"]
T.cal <- cal.cbs["1516", "T.cal"]
pnbd.ConditionalExpectedTransactions(params,
                                     T.star = 52,
                                     x,
                                     t.x,
                                     T.cal)
# [1] 25.45647
pnbd.PAlive(params,
            x,
            t.x,
            T.cal)
# [1] 0.997874

cet <- "pnbd.ConditionalExpectedTransactions"
for (i in seq(10, 25, 5)){
  cond.expectation <- match.fun(cet)(params,
                                     T.star = 52,
                                     x = i,
                                     t.x = 20,
                                     T.cal = 39)
  cat ("x:",i,"\t Expectation:",cond.expectation, fill = TRUE)
}
# x: 10 Expectation: 0.7062289
# x: 15 Expectation: 0.1442396
# x: 20 Expectation: 0.02250658
# x: 25 Expectation: 0.00309267

pnbd.PlotFrequencyInCalibration(params = params,
                                cal.cbs = cal.cbs,
                                censor = 7)

elog <- dc.SplitUpElogForRepeatTrans(elog)$repeat.trans.elog
x.star <- rep(0, nrow(cal.cbs))
cal.cbs <- cbind(cal.cbs, x.star)
elog.custs <- elog$cust
for (i in 1:nrow(cal.cbs)){
  current.cust <- rownames(cal.cbs)[i]
  tot.cust.trans <- length(which(elog.custs == current.cust))
  cal.trans <- cal.cbs[i, "x"]
  cal.cbs[i, "x.star"] <- tot.cust.trans - cal.trans
}
round(cal.cbs[1:3,], digits = 3)

T.star <- 39 # length of the holdout period
censor <- 7 # This censor serves the same purpose described above
x.star <- cal.cbs[,"x.star"]

comp <- pnbd.PlotFreqVsConditionalExpectedFrequency(params,
                                                    T.star,
                                                    cal.cbs,
                                                    x.star,
                                                    censor)

# pdf
# 2
rownames(comp) <- c("act", "exp", "bin")
round(comp, digits = 3)
# freq.0 freq.1 freq.2 freq.3 freq.4 freq.5 freq.6 freq.7+
# act 0.237 0.697 1.393 1.560 2.532 2.947 3.862 6.359
# exp 0.138 0.600 1.196 1.714 2.399 2.907 3.819 6.403
# bin 1411.000 439.000 214.000 100.000 62.000 38.000 29.000 64.000

