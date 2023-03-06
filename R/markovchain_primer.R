# https://spedygiorgio.github.io/markovchain/articles/an_introduction_to_markovchain_package.pdf

library(tidyverse)
library(markovchain)


# 3. The Markovchain package -----------------------------------------------


# 3.1 Creating markovchain objects -----------------------------------------



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



mcWeather <- new(Class            = "markovchain", 
                 states           = c("sunny", "cloudy", "rain"),
                 transitionMatrix = matrix(data = c(0.70, 0.2, 0.1,
                                                    0.3, 0.4, 0.3,
                                                    0.2, 0.45, 0.35), 
                                           byrow = byRow, 
                                           nrow = 3),
                 name             = "Weather")


defaultMc <- new("markovchain")

mcList <- new(Class        = "markovchainList", 
              markovchains = list(mcWeather, defaultMc),
              name         = "A list of Markov chains")


# 3.2 Handling markovchain objects -----------------------------------------

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
plot(mcWeather)


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

# 4. Probability with markovchain objects ---------------------------------


# 4.1 Conditional distributions -------------------------------------------

conditionalDistribution(mcWeather, "sunny")


# 4.2 Stationary states ---------------------------------------------------



steadyStates(mcWeather)

gamblerRuinMarkovChain <- function(moneyMax, prob = 0.5) {
  m <- matlab::zeros(moneyMax + 1)
  m[1,1] <- m[moneyMax + 1,moneyMax + 1] <- 1
  states <- as.character(0:moneyMax)
  rownames(m) <- colnames(m) <- states
  
  for(i in 2:moneyMax){
    m[i,i-1] <- 1 - prob
    m[i, i + 1] <- prob
  }
  
  new("markovchain", transitionMatrix = m,
      name = paste("Gambler ruin", moneyMax, "dim", sep = " "))
}

mcGR4 <- gamblerRuinMarkovChain(moneyMax = 4, prob = 0.5)
steadyStates(mcGR4)


#  4.3 Classification of states --------------------------------------------

absorbingStates(mcGR4)

absorbingStates(mcWeather)

P <- matlab::zeros(10)
P[1, c(1, 3)] <- 1/2;
P[2, 2] <- 1/3; P[2,7] <- 2/3;
P[3, 1] <- 1;
P[4, 5] <- 1;
P[5, c(4, 5, 9)] <- 1/3;
P[6, 6] <- 1;
P[7, 7] <- 1/4; P[7,9] <- 3/4;
P[8, c(3, 4, 8, 10)] <- 1/4;
P[9, 2] <- 1;
P[10, c(2, 5, 10)] <- 1/3;
rownames(P) <- letters[1:10]
colnames(P) <- letters[1:10]
probMc <- new("markovchain", transitionMatrix = P, name = "Probability MC")
summary(probMc)

transientStates(probMc)

probMcCanonic <- canonicForm(probMc)
probMc
probMcCanonic
is.accessible(object = probMc, from = "a", to = "c")
is.accessible(object = probMc, from = "g", to = "c")

E <- matrix(0, nrow = 4, ncol = 4)
E[1, 2] <- 1
E[2, 1] <- 1/3; E[2, 3] <- 2/3
E[3,2] <- 1/4; E[3, 4] <- 3/4
E[4, 3] <- 1

mcE <- new("markovchain", states = c("a", "b", "c", "d"),
           transitionMatrix = E,
           name = "E")
is.irreducible(mcE)

markovchain::period(mcE)

require(matlab)
mathematicaMatr <- zeros(5)
mathematicaMatr[1,] <- c(0, 1/3, 0, 2/3, 0)
mathematicaMatr[2,] <- c(1/2, 0, 0, 0, 1/2)
mathematicaMatr[3,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[4,] <- c(0, 0, 1/2, 1/2, 0)
mathematicaMatr[5,] <- c(0, 0, 0, 0, 1)
statesNames <- letters[1:5]
mathematicaMc <- new("markovchain", transitionMatrix = mathematicaMatr,
                      name = "Mathematica MC", states = statesNames)


# 4.4 First passage time distributions and means ---------------------------

# We conclude that the probability for the first rainy day to be the third one, 
# given that the current state is sunny, is given by:
firstPassagePdF <- firstPassage(object = mcWeather, state = "sunny", n = 10)
firstPassagePdF[3, "rain"]

firstPassagePdF[, "rain"] %>% plot(type = "l")

# To compute the mean first passage times, i.e. the expected number of days before it rains
# given that today is sunny, we can use the meanFirstPassageTime function:
meanFirstPassageTime(mcWeather)

# indicating e.g. that the average number of days of sun or cloud before rain is 
# 6.67 if we start counting from a sunny day, and 5 if we start from a cloudy day. 
# Note that we can also specify one or more destination states:
meanFirstPassageTime(mcWeather,"rain")

# We can check the result by averaging the first passage probability density function:
firstPassagePdF.long <- firstPassage(object = mcWeather, state = "sunny", n = 100)
sum(firstPassagePdF.long[,"rain"] * 1:100)


# 4.5 Mean recurrence time ------------------------------------------------

# The meanRecurrenceTime method gives the first mean recurrence time (expected number of
# steps to go back to a state if it was the initial one) for each recurrent state 
# in the transition probabilities matrix for a DTMC. Let’s see an example:
meanRecurrenceTime(mcWeather)

# 4.6 Absorption probabilities and mean absorption time -------------------


# 4.7 Committor probability -----------------------------------------------


# 4.8 Hitting probabilities -----------------------------------------------


# 5 Statistical analysis --------------------------------------------------



# 5.1 Simulation ----------------------------------------------------------

weathersOfDays <- rmarkovchain(n = 365, object = mcWeather, t0 = "sunny")
weathersOfDays[1:30]

patientStates <- rmarkovchain(n = 5, object = mcCCRC, t0 = "H", 
                              include.t0 = TRUE)
patientStates[1:10,]


# 5.2 Estimation ----------------------------------------------------------

createSequenceMatrix(stringchar = weathersOfDays,toRowProbs = TRUE)

weatherFittedMLE <- markovchainFit(data = weathersOfDays, method = "mle", name = "Weather")

weatherFittedMLE$estimate

weatherFittedMLE$standardError

weatherFittedMLE$confidenceLevel

weatherFittedMLE$lowerEndpointMatrix

weatherFittedMLE$upperEndpointMatrix

weatherFittedMLE$logLikelihood


# 5.3 Prediction ----------------------------------------------------------

# A 3-days forward predictions from markovchain object can be generated as follows, 
# assuming that the last two days were respectively “cloudy” and “sunny”.
predict(object = weatherFittedMLE$estimate, newdata = c("cloudy", "sunny"),
        n.ahead = 3)

# Given an initial two years health status, the 5-year ahead prediction of any CCRC guest is
predict(mcCCRC, newdata = c("H", "H"), n.ahead = 5)


# 5.4 Statistical Tests ---------------------------------------------------

# In this section, we describe the statistical tests: 
  # - ssessing the Markov property (verifyMarkovProperty),
  # - the order (assessOrder), 
  # - the stationary (assessStationarity) of a Markov chain sequence,
  # - and the divergence test for empirically estimated transition matrices (divergenceTest).

sample_sequence<-c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a",
                   "b", "a", "a", "b", "b", "b", "a")
verifyMarkovProperty(sample_sequence)

data(rain)
assessOrder(rain$rain)

assessStationarity(rain$rain, 10)

sequence <- c(0,1,2,2,1,0,0,0,0,0,0,1,2,2,2,1,0,0,1,0,0,0,0,0,0,1,1,
            2,0,0,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,2,1,0,
            0,2,1,0,0,0,0,0,0,1,1,1,2,2,0,0,2,1,1,1,1,2,1,1,1,1,1,1,1,1,1,0,2,
            0,1,1,0,0,0,1,2,2,0,0,0,0,0,0,2,2,2,1,1,1,1,0,1,1,1,1,0,0,2,1,1,
            0,0,0,0,0,2,2,1,1,1,1,1,2,1,2,0,0,0,1,2,2,2,0,0,0,1,1)
mc <- matrix(c(5/8, 1/4, 1/8, 1/4, 1/2, 1/4, 1/4, 3/8, 3/8), byrow = TRUE, nrow = 3)
rownames(mc) <- colnames(mc) <- 0:2

theoreticalMc <- as(mc, "markovchain")

verifyEmpiricalToTheoretical(data = sequence, object = theoreticalMc)

data(kullback)
verifyHomogeneity(inputList = kullback, verbose = TRUE)
