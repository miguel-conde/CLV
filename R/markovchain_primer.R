library(tidyverse)
library(markovchain)

data("weathersOfDays")


weatherStates <- c("sunny", "cloudy", "rain")
byRow <- TRUE
weatherMatrix <- matrix(data = c(0.70, 0.2, 0.1,
                                 0.3, 0.4, 0.3,
                                 0.2, 0.45, 0.35), 
                        byrow = byRow, 
                        nrow = 3,
                        dimnames = list(weatherStates, weatherStates))

mcWeather <- new(Class            = "markovchain", 
                 states           = weatherStates, 
                 byrow            = byRow,
                 transitionMatrix = weatherMatrix, 
                 name             = "Weather")



mcWeather <- new(class            = "markovchain", 
                 states           = c("sunny", "cloudy", "rain"),
                 transitionMatrix = matrix(data = c(0.70, 0.2, 0.1,
                                                    0.3, 0.4, 0.3,
                                                    0.2, 0.45, 0.35), 
                                           byrow = byRow, 
                                           nrow = 3),
                 name             = "Weather")


defaultMc <- new("markovchain")

mcList <- new(class        = "markovchainList", 
              markovchains = list(mcWeather, defaultMc),
              name         = "A list of Markov chains")

initialState <- c(0, 1, 0)
after2Days   <- initialState * (mcWeather * mcWeather)
after7Days   <- initialState * (mcWeather ^ 7)
after2Days
after7Days


fvals <- function(mchain, initialstate, n) {
  out <- data.frame()
  names(initialstate) <- names(mchain)
  for (i in 0:n)
  {
    iteration <- initialstate*mchain^(i)
    out <- rbind(out, iteration)
  }
  out <- cbind(out, i = seq(0, n))
  out <- out[, c(4, 1:3)]
  
  return(out)
}

fvals(mchain = mcWeather, initialstate = c(90 ,5, 5), n = 4)

states(mcWeather)
names(mcWeather)
dim(mcWeather)
name(mcWeather)

markovchain:::sort(mcWeather)
transitionProbability(mcWeather, "cloudy", "rain")
mcWeather[2,3]
mcWeather["cloudy", "rain"]
print(mcWeather)
show(mcWeather)


mcDf <- as(mcWeather, "data.frame")
mcDf
mcNew <- as(mcDf, "markovchain")
mcNew
name(mcNew) <- "new weather"
mcNew

mcIgraph <- as(mcWeather, "igraph")
if (requireNamespace("msm", quietly = TRUE)) {
  require(msm)
  Q <- rbind ( c(0, 0.25, 0, 0.25),
               c(0.166, 0, 0.166, 0.166),
               c(0, 0.25, 0, 0.25),
               c(0, 0, 0, 0) )
  cavmsm <- msm(state ~ years, subject = PTNUM, data = cav, qmatrix = Q, death = 4)
  msmMc <- as(cavmsm, "markovchain")
  msmMc
} else {
  message("msm unavailable")
}

if (requireNamespace("etm", quietly = TRUE)) {
  library(etm)
  data(sir.cont)
  sir.cont <- sir.cont[order(sir.cont$id, sir.cont$time), ]
  for (i in 2:nrow(sir.cont)) {
    if (sir.cont$id[i]==sir.cont$id[i-1]) {
      if (sir.cont$time[i]==sir.cont$time[i-1]) {
        sir.cont$time[i-1] <- sir.cont$time[i-1] - 0.5
      }
    }
  }
  tra <- matrix(ncol=3,nrow=3,FALSE)
  tra[1, 2:3] <- TRUE
  tra[2, c(1, 3)] <- TRUE
  tr.prob <- etm::etm(sir.cont, c("0", "1", "2"), tra, "cens", 1)
  tr.prob
  etm2mc<-as(tr.prob, "markovchain")
  etm2mc
} else {
  message("etm unavailable")
}


myMatr<-matrix(c(.1,.8,.1,.2,.6,.2,.3,.4,.3), byrow=TRUE, ncol=3)
myMc<-as(myMatr, "markovchain")
myMc


stateNames = c("H", "I", "D")
Q0 <- new("markovchain", states = stateNames,
          transitionMatrix =matrix(c(0.7, 0.2, 0.1,0.1, 0.6, 0.3,0, 0, 1),
                                   byrow = TRUE, nrow = 3), name = "state t0")
Q1 <- new("markovchain", states = stateNames,
          transitionMatrix = matrix(c(0.5, 0.3, 0.2,0, 0.4, 0.6,0, 0, 1),
                                    byrow = TRUE, nrow = 3), name = "state t1")
Q2 <- new("markovchain", states = stateNames,
          transitionMatrix = matrix(c(0.3, 0.2, 0.5,0, 0.2, 0.8,0, 0, 1),
                                    byrow = TRUE,nrow = 3), name = "state t2")
Q3 <- new("markovchain", states = stateNames,
          transitionMatrix = matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 1),
                                    byrow = TRUE, nrow = 3), name = "state t3")
mcCCRC <- new("markovchainList",markovchains = list(Q0,Q1,Q2,Q3),
              name = "Continuous Care Health Community")
print(mcCCRC)


mcCCRC[[1]]
dim(mcCCRC)

# Probability with markovchain objects ------------------------------------

conditionalDistribution(mcWeather, "sunny")
