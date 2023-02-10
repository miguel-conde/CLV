library(tidyverse)

# Tasa de descuento
d <- .10

# Presente
n_actual <- 37

# Pasado  
n_past <- seq(1, n_actual-1)

# Futuro
m <- n_actual + 24
n_future <- seq(n_actual+1, m)

ns <- c(n_past, n_actual, n_future)

disc_past_present <- (1-d)^((n_actual-1):0)

plot(c(n_past, n_actual), disc_past_present)


disc_future <- (1-d)^-(n_future - n_actual)
plot(n_future, disc_future)


plot(c(n_past, n_actual, n_future), c(disc_past_present, disc_future), type = "o")
abline(v = n_actual, lty = 2)
abline(h = 1, lty = 2)
