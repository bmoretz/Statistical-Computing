library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(GGally)
library(ggthemes)
library(scales)
library(reshape2)
library(skimr)
library(gridExtra)

#####################################################################
######################### Chapter 2 Lab #############################
#####################################################################

theme_set(theme_sjplot())

path.data <- "D:/Projects/MSDS-RiskAnalytics/datasets"
setwd(path.data)

# 2.4.1
# Data Analytics

dat <- read.csv("Stock_bond.csv", header = T)

names(dat)
attach(dat)

par(mfrow = c(1, 2))

plot(GM_AC, type = "l")
plot(F_AC, type = "l")

n <- dim(dat)[1]
GMReturn <- GM_AC[-1] / GM_AC[-n] - 1
FReturn <- F_AC[-1] / F_AC[-n] - 1

par(mfrow = c(1, 1))
plot(GMReturn, FReturn)

# Problems

# 1.)
# Do GM and Ford returns seem positively correlated? Do you notice any outliying returns? 

cor(GMReturn, FReturn) # .61

# If "yes," do outlying GM returns seem to occur with outliying Ford returns?

logGMReturn <- diff(log(GM_AC))

plot(GMReturn, logGMReturn)
cor(GMReturn, logGMReturn)

# MSFT / MRK

MSReturn <- MSFT_AC[-1] / MSFT_AC[-n] - 1
logMSReturn <- diff(log(MSFT_AC))

MRKReturn <- MRK_AC[-1] / MRK_AC[-n] - 1
logMRKReturn <- diff(log(MRK_AC))

plot(MSReturn, MRKReturn)

cor(MSReturn, MRKReturn)

plot(MSReturn, logMSReturn)

# 2.4.2
# Simulations

niter <- 1e5 # number of iterations
set.seed(2009) # reproducible

seed.capital <- 5e4
initial.investment <- log(1e6)
profit.threshold <- log(1.1e6)
loss.threshold <- log(9.5e5)

target.profit <- 1e5

ret.avg <- 0.05
ret.sd <- 0.23
market.open <- 253

simulate_market <- function(days) {
  # generate random returns for N days
  r <- rnorm(days, mean = ret.avg / market.open,
            sd = ret.sd / sqrt(market.open))

  cumsum(r) # return the final log price after N days.
}

# setup storage
outcomes <- list(below = rep(0, niter))

# Simulation: Probability dips below $950,000.
for (i in 1:niter) {

  logPrice = initial.investment + simulate_market(45) # simulate 45 trading days.

  minlogP = min(logPrice) # miniumum price over next 45 days

  outcomes$below[i] = as.numeric(minlogP < loss.threshold)
}

print(paste0("Probability the value of the stock is below $950,000 at least one of next 45 sessions: ", round(mean(outcomes$below), 3) * 100, "%"))

# reset seed for next simulation.

# Suppose the hedge fund will sell the stock for a profit of at least $100,000 if the value of the stock rises to at least 
# $1,100,000 at the end of one of the first 100 trading days, sell it for a sloss if the value falls below $950,000 at the end of
# one of the first 100 trading days, or sell it (for "FMV") after 100 trading days if the closing price has stayed between $950,000 and $1,100,000.

set.seed(2009) # reproducible

outcomes <- list(above = rep(0, niter),
                 below = rep(0, niter),
                 middle = rep(0, niter),
                 pnl = rep(0, niter),
                 ret = rep(0, niter))

for (i in 1:niter) {

  logPrice = initial.investment + simulate_market(100) # simulate 100 trading days.

  suppressWarnings({
    # ignore Inf returned if condition not meet.
    profit.day <- min(which(logPrice >= profit.threshold))
    loss.day <- min(which(logPrice <= loss.threshold))
  })

  is.market <- profit.day == Inf && loss.day == Inf

  # What was the exit condition of the position, hince the final price of the stock?
  daysOpen <- ifelse(is.market, length(logPrice),
                       min(profit.day, loss.day))

  outcomes$above[i] <- min(profit.day) < min(loss.day)
  outcomes$middle[i] <- is.market
  outcomes$below[i] <- min(loss.day) < min(profit.day)

  # p&l = ending value - initial investment
  pnl <- exp(logPrice[daysOpen]) - exp(initial.investment)

  outcomes$pnl[i] <- ifelse(is.market, pnl,
                            ifelse(pnl >= 0, target.profit, -seed.capital))

  outcomes$ret[i] <- (outcomes$pnl[i] / seed.capital) / daysOpen
}

stopifnot(sum(outcomes$above) + sum(outcomes$below) + sum(outcomes$middle) == niter)

prob.profit.target <- sum(outcomes$pnl >= target.profit) / length(outcomes$pnl) # Probability of profit over $100,000.
print(paste0("Probability the hedge fund (strategy) returns over $100,000 in profit: ", round(mean(prob.profit.target), 3) * 100, "%"))

prob.loss <- sum(outcomes$below) / length(outcomes$below) # Probability of loss
print(paste0("Probability the hedge fund (strategy) returns a loss: ", round(mean(prob.loss), 3) * 100, "%"))

ev.pnl <- mean(outcomes$pnl) # Expected P&L
print(paste0("Expected profit/loss of the hedge fund (strategy): $", round(ev.pnl, 2)))

ev.ret <- mean(outcomes$ret) # Expected Return (TW)
print(paste0("Expected (time-weighted) return of the hedge fund (strategy): ", round(ev.ret, 5) * 100, "%"))

# 2.4.3
# Simulating Geometric Random Walk

set.seed(2012)
n = 253
par(mfrow = c(3, 3))
for (i in (1:9)) {
  logr = rnorm(n, 0.05 / 253, 0.2 / sqrt(253))
  price = c(120, 120 * exp(cumsum(logr)))
  plot(price, type = "b")
}

# What are the mean & sd of the log-returns for 1 year?

mean(logr) * 100
sd(logr)

# 2.4.4
# McDonald's Stock

par(mfrow = c(1, 1))

data <- read.csv("MCD_PriceDaily.csv", header = T)
head(data)
adjPrice <- data[, 7]

n <- length(adjPrice)

rets <- adjPrice[2:n] / adjPrice[1:(n - 1)] - 1
lrets <- log(adjPrice[2:n] / adjPrice[1:(n - 1)])

plot(rets, lrets, type = 'p', xlab = 'return', ylab = 'log return')
grid()
print(cor(rets, lrets))

print(sprintf('returns: mean= %f; std= %f', mean(rets), sd(rets)))
print(sprintf('log-returns: mean= %f; std= %f', mean(lrets), sd(lrets)))

t.test(x = rets, y = lrets, paired = TRUE)

