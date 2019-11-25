
# Suppose that the daily log returns on a stock are i.i.d. and normally distributed with mean 0.001 and a sd of 0.015.
# Suppose you buy $1,000 worth of this stock.

# What is the probability that after one trading day your investment is worth less than $990?

days <- 1
ret.avg <- 0.001
ret.sd <- 0.015

price.threshold <- log(990 / 1000)

round(pnorm(price.threshold, mean = ret.avg * days, sd = ret.sd * sqrt(days)), 4) * 100

# What is the probability after 5 trading days your investment is worth less than $990?

days <- 5

round(pnorm(price.threshold, mean = ret.avg * days, sd = ret.sd * sqrt(days)), 4) * 100

# The yearly log returns on a stock are normally distributed with mean 0.1 and sd 0.2.

ret.avg <- 0.1
ret.sd <- 0.2

# The stock is selling at $100 today.

# What is the probability that 1 year from now it is selling at $110 or more?

days <- 252
price.threshold <- log(110 / 100)

round(pnorm(price.threshold, mean = ret.avg, sd = ret.sd, lower.tail = F), 5) * 100

# The yearly log returns on a stock are normally distributed with mean 0.08 and sd 0.15.

years <- 2

ret.avg <- 0.08
ret.sd <- 0.15

# What is the probability that 2 years from now it is selling at $90 or more?
price.threshold <- log(90 / 80)

round(pnorm(price.threshold, mean = ret.avg * years, sd = ret.sd * sqrt(years), lower.tail = F), 5) * 100

# Suppose the prices of a stock at times 1, 2 and 3 are:

prices <- c(95, 103, 98)

# Find R(3):

ret <- (prices[3] - prices[1]) / prices[1]
log(1 + ret)

