######################################################
# Module 1, Monte Carlo Simulations
# Moretz, Brandon
######################################################

# Suppose a hedge fund owns $1,000,000 of stock and used $50,000 of its own capital and $950,000 in borrowed money for the purchase.
# Suppose that if the value of the stock falls below $950,000 at the end of any trading day, then the hedge fund will sell all the stock and repay the loan.
# This will wipe out its $50,000 investment.
# The hedge fund is said to be leveraged 20:1 since its position is 220 times the amount of its own capital invested.

# Suppose that the daily log returns on the stock have a mean of 0.05/year and a standard deviation of 0.23/year. 
# These can be converted to rates per trading day by dividing by 253 and sqrt(253), respectively.

# What is the probability that the value of the stock will be below $950,000 at the close of at least one of the next 45 trading days?

niter <- 1e5 # number of iterations
set.seed(2009) # reproducible

seed.capital <- 5e4
initial.investment <- log(1e6)
profit.threshold <- log(1.1e6)
loss.threshold <- log(9.5e5)
target.profit <- 1e5

simulate_market <- function(days) {
  # generate random returns for N days
  r <- rnorm(days, mean = 0.05 / 253,
             sd = 0.23 / sqrt(253))
  
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
# reproducible
set.seed(2009)

outcomes <- list(above = rep(0, niter),
                 below = rep(0, niter),
                 middle = rep(0, niter),
                 pnl = rep(0, niter),
                 capped_pnl = rep(0, niter),
                 ret = rep(0, niter))

for (i in 1:niter) {
  
  # simulate 100 trading days.
  logPrice = initial.investment + simulate_market(100)
  
  suppressWarnings({
    # ignore Inf returned if condition not meet.
    profit.day <- min(which(logPrice >= profit.threshold))
    loss.day <- min(which(logPrice <= loss.threshold))
  })
  
  is.market <- profit.day == Inf && loss.day == Inf
  
  # What was the exit condition of the position, hince the final price of the stock?
  days.open <- ifelse(is.market, length(logPrice),
                      min(profit.day, loss.day))
  
  outcomes$above[i] <- min(profit.day) < min(loss.day)
  outcomes$middle[i] <- is.market
  outcomes$below[i] <- min(loss.day) < min(profit.day)
  
  # p&l = ending value - initial investment
  pnl <- exp(logPrice[days.open]) - exp(initial.investment)
  
  # market pnl = use FMV, otherwise cap p/l
  outcomes$pnl[i] <- pnl
  
  outcomes$capped_pnl[i] <- ifelse(is.market, pnl, 
                                   ifelse(pnl > 0, target.profit, -seed.capital))
  
  # Calculate return (time-weighted)
  outcomes$ret[i] <- (outcomes$pnl[i] / seed.capital) / days.open
}

# Verify we captured every simulation outcome.
stopifnot(sum(outcomes$above) + sum(outcomes$below) + sum(outcomes$middle) == niter)

mean(outcomes$pnl) # 8910.54
mean(outcomes$capped_pnl) # 9922.63

prob.profit <- sum(outcomes$above) / length(outcomes$pnl)
print(paste0("Probability the hedge fund (strategy) returns over $100,000 in profit: ", round(mean(prob.profit), 3) * 100, "%"))

prob.loss <- sum(outcomes$below) / length(outcomes$below)
print(paste0("Probability the hedge fund (strategy) returns a loss: ", round(mean(prob.loss), 3) * 100, "%"))

ev.pnl <- mean(outcomes$pnl)
print(paste0("Expected profit/loss of the hedge fund (strategy): $", round(ev.pnl, 2)))

ev.ret <- mean(outcomes$ret)
print(paste0("Expected (time-weighted) return of the hedge fund (strategy): ", round(ev.ret, 5) * 100, "%"))
