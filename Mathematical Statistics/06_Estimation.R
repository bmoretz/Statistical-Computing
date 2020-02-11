# Suppose X1, X2, ... Xn are a random sample from the Cauchy distribution with pdf f(x;theta) = 1/(pi(1 + (x - theta)^2))

x <- c(1, 2, 2, 3)

g <- function(theta) sum(log( 1 + (x-theta)^2))

obj <- optimize(g, interval = c(0, 4))

obj$minimum

obj$objective

logL <- function(theta) sum(log(dcauchy(x, theta)))

optimize(logL, interval = c(0, 4), maximum = T)
