# https://www.analyticsvidhya.com/blog/2020/10/a-definitive-guide-for-predicting-customer-lifetime-value-clv/

library(tidyverse)
library(readxl)
library(lubridate)

source(here::here("R", "global.R"))

PROFIT_MARGIN <- 0.05

# DATA --------------------------------------------------------------------

# https://archive.ics.uci.edu/ml/datasets/online+retail
#
# InvoiceNo: Invoice number. Nominal, a 6-digit integral number uniquely assigned to each transaction. If this code starts with letter 'c', it indicates a cancellation.
# StockCode: Product (item) code. Nominal, a 5-digit integral number uniquely assigned to each distinct product.
# Description: Product (item) name. Nominal.
# Quantity: The quantities of each product (item) per transaction. Numeric.
# InvoiceDate: Invice Date and time. Numeric, the day and time when each transaction was generated.
# UnitPrice: Unit price. Numeric, Product price per unit in sterling.
# CustomerID: Customer number. Nominal, a 5-digit integral number uniquely assigned to each customer.
# Country: Country name. Nominal, the name of the country where each customer resides.

online_retail_raw <- read_excel(FILE_CLTV_DATA) 

online_retail <- online_retail_raw %>% 
  janitor::clean_names() %>% 
  mutate(customer_id = as.integer(customer_id),
         total_sales = unit_price * quantity) %>% 
  select(customer_id, invoice_no, invoice_date, quantity, unit_price, total_sales) %>% 
  drop_na() %>% 
  filter(!str_detect(invoice_no, "^C"))

online_retail

summary(online_retail)

range(online_retail$invoice_date)
online_retail$customer_id %>% unique() %>% length()
online_retail$quantity %>% sum()
online_retail$total_sales %>% sum()

# HISTORICAL APPROACH -----------------------------------------------------


# Aggregate Model ---------------------------------------------------------

customer <- online_retail %>% 
  mutate(invoice_date = as_date(invoice_date)) %>% 
  group_by(customer_id) %>% 
  summarise(age = as.numeric(max(invoice_date) - min(invoice_date)),
            frequency = n(),
            total_sales = sum(total_sales))

calc_cltv <- function(x) {
  avg_sales      <- mean(x$total_sales)
  purchase_freq  <- mean(x$frequency)
  retention_rate <- nrow(x %>% filter(frequency > 1)) / nrow(x)
  churn          <- 1 - retention_rate
  
  CLTV <- avg_sales * purchase_freq / churn * PROFIT_MARGIN
  return(CLTV)
}

CLTV <- calc_cltv(customer)
CLTV

customer$total_sales %>% summary()
customer$total_sales %>% hist(breaks = 200)

# Cohort Model ------------------------------------------------------------

# So, the question would be- “how do we group the customers?”

# The most common way to group customers into cohorts is by the start date of a 
# customer, typically by month. The best choice will depend on the customer 
# acquisition rate, seasonality of the business, and whether additional customer 
# information can be used.

customer <- online_retail %>% 
  mutate(invoice_date = as_date(invoice_date)) %>% 
  group_by(customer_id) %>% 
  summarise(start_month = month(min(invoice_date)),
            frequency = n(),
            total_sales = sum(total_sales))

CLTV <- customer %>% 
  group_by(start_month) %>% 
  group_map(~ calc_cltv(.x), .keep = TRUE)
names(CLTV) <- month.name
CLTV <- CLTV %>% bind_rows() %>% gather(month, cltv)
CLTV

# PREDICTIVE APPROACH -----------------------------------------------------


# Machine Learning Model --------------------------------------------------


# Probabilistic Model -----------------------------------------------------

# https://www.clvtools.com/articles/CLVTools.html
library("CLVTools")

customer <- online_retail %>% 
  mutate(invoice_date = as_date(invoice_date)) %>% 
  group_by(customer_id, invoice_date) %>% 
  summarise(total_sales = sum(total_sales), .groups = "drop")

clv_customer <- clvdata(customer,  
                        date.format="ymd", 
                        time.unit = "day",
                        # estimation.split = 40,
                        name.id = "customer_id",
                        name.date = "invoice_date",
                        name.price = "total_sales")
clv_customer

summary(clv_customer)

plot(clv_customer, which = "tracking") # Default
plot(clv_customer, which = "frequency")
plot(clv_customer, which = "spending")
plot(clv_customer, which = "interpurchasetime")

est_bgnbd <- bgnbd(clv.data = clv_customer)
est_pnbd <- pnbd(clv.data = clv_customer)

est_bgnbd

plot(est_bgnbd, which = "tracking") # Default)
plot(est_bgnbd, which = "pmf")

summary(est_bgnbd)
coef(est_bgnbd)
confint(est_bgnbd)
logLik(est_bgnbd)
vcov(est_bgnbd)

results <- predict(est_bgnbd, prediction.end = 10)
results


est_gg<- gg(clv.data = clv_customer)

est_gg

summary(est_gg)

plot(est_gg)

results_spending <- predict(est_gg)

results_spending


###
predict(est_bgnbd, prediction.end = 10)
predict(est_bgnbd, prediction.end = 10, predict.spending = TRUE)
predict(est_bgnbd, prediction.end = 10, predict.spending = est_gg)

predict(est_pnbd, prediction.end = 10)
predict(est_pnbd, prediction.end = 10, predict.spending = TRUE)
predict(est_pnbd, prediction.end = 10, predict.spending = est_gg)
