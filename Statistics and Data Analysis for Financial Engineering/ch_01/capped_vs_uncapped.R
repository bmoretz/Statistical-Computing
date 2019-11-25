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

# reset seed for next simulation.

# Suppose the hedge fund will sell the stock for a profit of at least $100,000 if the value of the stock rises to at least 
# $1,100,000 at the end of one of the first 100 trading days, sell it for a sloss if the value falls below $950,000 at the end of
# one of the first 100 trading days, or sell it (for "FMV") after 100 trading days if the closing price has stayed between $950,000 and $1,100,000.
# reproducible
set.seed(2009)

outcomes <- list(above = rep(0, niter),
                 below = rep(0, niter),
                 middle = rep(0, niter),
                 uncapped_pnl = rep(0, niter),
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
  outcomes$uncapped_pnl[i] <- pnl
  
  outcomes$capped_pnl[i] <- ifelse(is.market, pnl, 
                                   ifelse(pnl > 0, target.profit, -seed.capital))
  
  # Calculate return (time-weighted)
  outcomes$ret[i] <- (outcomes$capped_pnl[i] / seed.capital) / days.open
}

# Verify we captured every simulation outcome.
stopifnot(sum(outcomes$above) + sum(outcomes$below) + sum(outcomes$middle) == niter)

mean(outcomes$uncapped_pnl) # 8910.544
mean(outcomes$capped_pnl) # 9922.63
