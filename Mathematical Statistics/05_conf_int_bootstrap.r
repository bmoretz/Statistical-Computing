library(data.table)
library(ggplot2)
library(here)
library(boot) # alternate to verify long-hand calcs

data.dir <- paste0(here::here(), "/datasets/")

# gamma example

my.sample <- rgamma(16, 1, 1/2)
N <- 10e5

my.boot <- numeric(N)

for(i in 1:N)
{
  x <- sample(x = my.sample, size = length(my.sample), replace = T) 
  my.boot[i] <- mean(x)
}

ggplot(data.table(value = my.sample), aes(value, fill = ..count..)) +
  geom_histogram(bins = 8)

ggplot(data.table(value = my.boot), aes(value, fill = ..count..)) +
  geom_histogram(bins = 30)

mean(my.boot)
sd(my.boot)

### Bangladesh

Bangladesh <- data.table(read.csv(paste0(data.dir, "Bangladesh.csv"),
                                  header = T))

Arsenic <- Bangladesh$Arsenic

mean(Arsenic) #  sample mean

ggplot(data.table(value = Arsenic)) +
  geom_histogram(aes(value, fill = ..count..), bins = 30) +
  scale_x_continuous(labels = comma) +
  labs(title = "Arsenic")

ggplot(data.table(value = Arsenic), aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Arsenic")

n <- length(Arsenic)
N <- 10e3

arsenic.mean <- numeric(N)

for(i in 1:N)
{
  x <- sample(Arsenic, n, replace = T)
  arsenic.mean[i] <- mean(x)
}

ggplot(data.table(value = arsenic.mean)) +
  geom_histogram(aes(value, fill = ..count..), bins = 30) +
  geom_vline(xintercept = mean(arsenic.mean), col = "darkorange", alpha = .6, lwd = 2) +
  scale_y_continuous(labels = comma) +
  labs(title = "Arsenic Bootstrap Distribution of Means")

ggplot(data.table(value = arsenic.mean), aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Arsenic")

mean(arsenic.mean) # Bootstrap mean
mean(arsenic.mean) - mean(Arsenic) # bias
sd(arsenic.mean) # bootstrap SE

# this is not accurate: CLT cannot be used here
lq <- mean(arsenic.mean) - 1.96 * sd(arsenic.mean)
uq <- mean(arsenic.mean) + 1.96 * sd(arsenic.mean)

sum(arsenic.mean < lq) / N
sum(arsenic.mean > uq) / N

alpha <- .05
quantile(arsenic.mean, c(alpha/2, 1 - alpha/2)) # confidence intervals

boot.fn <- function(data, index) {
  mean(data[index])
}

arsenic.boot <- boot(Arsenic, boot.fn, R = N)

### NC birth weights

NCBirths <- data.table(read.csv(paste0(data.dir, "NCBirths2004.csv"),
                                  header = T))

weights <- NCBirths$Weight

boot.fn <- function(data, index) {
  mean(data[index])  
}

bw.boot <- boot(weights, boot.fn, R = 1000)

# long-hand

N <- 1e5
bw.mean <- numeric(N)

for(i in 1:N)
{
  x <- sample(weights, size = length(weights), replace = T)
  bw.mean[i] <- mean(x)
}

alpha <- 0.05
conf <- quantile(bw.mean, c(Lower = alpha/2, Upper = 1 - alpha/2))

# NC Birth weights, bootstrapped

ggplot(data.table(values = bw.mean), aes(values)) +
  geom_histogram(aes(fill = ..count..), bins = 30) +
  geom_vline(xintercept = mean(weights), col = "darkorange", alpha = .4, lwd = 1) +
  geom_vline(xintercept = conf[1], col = "cornflowerblue") +
  geom_vline(xintercept = conf[2], col = "cornflowerblue") +
  geom_rug(col = "darkred") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)

### testosterone study

Skateboard <- data.table(read.csv(paste0(data.dir, "Skateboard.csv"),
                                  header = T))

testF <- Skateboard[Experimenter == "Female"]$Testosterone
testM <- Skateboard[Experimenter == "Male"]$Testosterone

nf <- length(testF)
nm <- length(testM)

N <- 10e4
TestMean <- numeric(N)

for(i in 1:N)
{
  sampleF <- sample(testF, nf, replace = T)
  sampleM <- sample(testM, nm, replace = T)
  
  TestMean[i] <- mean(sampleF) - mean(sampleM)
}

ggplot(data.table(result = TestMean), aes(result)) +
  geom_histogram(aes(fill = ..count..), bins = 30) +
  geom_vline(xintercept = mean(testF) - mean(testM), col = "darkorange", lty = 2, lwd = 1.5, alpha = .4) +
  labs(title = "Bootstrap distribution of difference in means")

ggplot(data.table(value = TestMean), aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Testosterone")

# bootstrap statistics

mean(testF) - mean(testM)
mean(TestMean)

sd(TestMean)

quantile(TestMean, c(0.025, .975)) # confidence intervals

mean(TestMean) - ( mean(testF - mean(testM)) ) # bias

observed <- mean(testF) - mean(testM)

result <- numeric(N)

for(i in 1:N)
{
  index <- sample(nf + nm, nf, replace = F)
  result[i] <- mean(Skateboard[index]$Testosterone) - mean(Skateboard[-index]$Testosterone)
}

p <- (sum(result >= observed) + 1) / (N + 1)

ggplot(data.table(result), aes(result)) +
  geom_histogram(aes(fill = ..count..), bins = 30) +
  geom_vline(xintercept = observed, col = "darkorange", lty = 2, lwd = 1.5, alpha = .4) +
  labs(title = "Permutation Test: Difference in means")

ggplot(data.table(value = result), aes(sample = value)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Testosterone Permutation Test")

### Verizon data

# Is the average repair time lower for ILEC lower than CLEC customers?

Verizon <- data.table(read.csv(paste0(data.dir, "Verizon.csv"),
                               header = T))

Time.ILEC <- Verizon[Group == "ILEC"]$Time
Time.CLEC <- Verizon[Group == "CLEC"]$Time

observed <- mean(Time.ILEC) - mean(Time.CLEC)

ggplot(data.table(values = Time.ILEC), aes(values)) +
  geom_histogram(aes(y = ..density.., fill = ..count..), bins = 30) +
  geom_density(aes(y = ..density..), col = "darkorange")

ggplot(data.table(values = Time.CLEC), aes(values)) +
  geom_histogram(aes(y = ..density.., fill = ..count..), bins = 30) +
  geom_density(aes(y = ..density..), col = "darkorange")

# bootstrap difference in means

set.seed(123)

N <- 10e4
results <- numeric(N)

N.ILEC <- length(Time.ILEC)
N.CLEC <- length(Time.CLEC)

for(i in 1:N)
{
  ilec <- sample(Time.ILEC, N.ILEC, replace = T)
  clec <- sample(Time.CLEC, N.CLEC, replace = T)
  
  results[i] <- mean(ilec) - mean(clec)
}

alpha <- 0.05
observed - mean(results) # Bias

sd(results)

quantile(results, c(Lower = alpha/2, Upper = 1 - alpha/2))

aov(Time ~ Group, data = Verizon)

# ratio of means test

observed <- mean(Time.ILEC) / mean(Time.CLEC)

N <- 10e3
time.ratio.mean <- numeric(N)
for(i in 1:N)
{
  ILEC.sample <- sample(Time.ILEC, length(Time.ILEC), replace = T)
  CLEC.sample <- sample(Time.CLEC, length(Time.CLEC), replace = T)
  
  time.ratio.mean[i] <- mean(ILEC.sample) / mean(CLEC.sample)
}

ggplot(data.table(result), aes(result)) +
  geom_histogram(aes(fill = ..count..), bins = 30) +
  geom_vline(xintercept = observed, col = "darkorange", lty = 2, lwd = 1.5, alpha = .4) +
  labs(title = "Permutation Test: Difference in means")

result <- numeric(N)

for(i in 1:N)
{
  index <- sample(nrow(Verizon), length(Time.CLEC), replace = F)
  result[i] <- mean(Verizon[index]$Time) / mean(Verizon[-index]$Time)
}

observed - mean(time.ratio.mean) # bias
sd(time.ratio.mean) # standard error

quantile(time.ratio.mean, c(0.025, .975))

time.ratio.bias <- mean(time.ratio.mean) - mean(Time.ILEC)/mean(Time.CLEC) # Bias

time.ratio.bias / sd(time.ratio.mean)


# Significance
p <- (sum(result >= observed) + 1) / ( N + 1)

ggplot(data.table(result), aes(result)) +
  geom_histogram(aes(fill = ..count..), bins = 30) +
  geom_vline(xintercept = observed, col = "darkorange", lty = 2, lwd = 1.5, alpha = .4) +
  scale_y_continuous(labels = comma) +
  labs(title = "Permutation Test: Difference in means")

