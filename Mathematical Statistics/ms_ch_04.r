N <- 1000
Xbar <- numeric(N) # space for results

for( i in 1:n)
{
  x <- rexp(100, rate = 1/15) # draw random sample of size 100
  Xbar[i] <- mean(x)
}

hist(Xbar)

qqnorm(Xbar)
qqline(Xbar)

mean(Xbar)
sd(Xbar)

maxY <- numeric(N)

for( i in 1:N)
{
  y <- runif(12) # draw random sample of size 12
  maxY[i] <- max(y)
}

hist(maxY)

X <- rpois(10^4, 5) # Draw 10^4 values from Pois(5)
Y <- rpois(10^4, 12) # Draw 10^4 values from Pois(12)

W <- X + Y

hist(W, prob = T) # prob = T, scales hist to 1
lines(2:35, dpois(2:35, 17), type = "b")

mean(W)
var(W)

sd(W) / sqrt(N)

Xbar <- numeric(N)

for( i in 1:N )
{
  x <- rgamma(30, shape = 5, rate = 2)
  Xbar[i] <- mean(x)
}

hist(Xbar)

qqnorm(Xbar)
qqline(Xbar)

mean(Xbar)
sd(Xbar)

mean(Xbar > 3)

# work

p <- c(3, 4, 6, 6)

sample.size <- 2
perm <- permutations( length(p), sample.size, repeats.allowed = T)

d <- data.table( c1 = p[perm[, 1]], c2 = p[perm[, 2]])
d[, m := (c1 + c2) / 2]

( 3 - 5/2 ) / (sqrt(5/2^2) / sqrt(30))

pnorm(2.44949, lower.tail = F)


z <- (0.53333 - 0.5) / 0.0289
pbinom(z, prob = .5, size = 300, lower.tail = T)

pbinom(160, size = 300, prob = .5, lower.tail = T)


