

D:\Projects\Statistical-Computing

qnorm(.5)


x <- seq(-1, 1, 0.01)
y <- qnorm(x)

plot(x, y)

plot(x, dnorm(x))

?pnorm

f <- function(x, mu = 0, sgma = 1) (1 / sgma*sqrt(2*pi)) ^ exp( (-1/2) * ( (x - mu)/sgma)**2 )

par(mfrow=c(1,1))
plot(x, f(x))
segments(a, 0, a, f(a), col = "blue")

qnorm(.455)
pnorm(-.11)


x <- c(17.7, 22.6, 26.1, 28.3, 30, 31.2, 31.5, 33.5, 34.7, 36)

qqnorm(x)
qqline(x)

x <- c(3, 6, 15, 15, 17, 19, 24)
plot.ecdf(x)

x <- rnorm(25)
plot.ecdf(x, xlim = c(-4, 4))
curve(pnorm(x), col = "blue", add = T)

Beerwings <- read.csv("http://sites.google.com/site/chiharahesterberg/data2/Beerwings.csv")

plot(Beer ~ Hotwings, data = Beerwings, col = Gender)

## Ex.1

x <- c(3, 5, 8, 15, 20, 21, 24)

mean(x)
median(x)

lnx <- log(x)

mean(lnx)
median(lnx)

# Ex. 2
x <- c(1, 2, 4, 5, 6, 8, 11, 15)

mean(x)
median(x)

xbar <- sqrt(x)

mean(xbar)
median(xbar)

Flights <- 