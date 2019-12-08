library(data.table)
library(ggplot2)

data.dir <- "D:/Projects/Statistical-Computing/datasets/"

### Bangladesh

Bangladesh <- data.table(read.csv(paste0(data.dir, "Bangladesh.csv"),
                                 header = T))

Arsenic <- Bangladesh$Arsenic

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

uq <- mean(arsenic.mean) - 1.96 * sd(arsenic.mean)
lq <- mean(arsenic.mean) + 1.96 * sd(arsenic.mean)

sum(arsenic.mean > lq) / N
sum(arsenic.mean < uq) / N

quantile(arsenic.mean, c(0.025, 0.975)) # confidence intervals

### testosterone study

Skateboard <- data.table(read.csv(paste0(data.dir, "Skateboard.csv"),
                                  header = T))

testF <- Skateboard[Experimenter == "Female"]$Testosterone
testM <- Skateboard[Experimenter == "Male"]$Testosterone

nf <- length(testF)
nm <- length(testM)

N <- 10e3
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

# Statistics

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

Verizon <- data.table(read.csv(paste0(data.dir, "Verizon.csv"),
                               header = T))

Time.ILEC <- Verizon[Group == "ILEC"]$Time
Time.CLEC <- Verizon[Group == "CLEC"]$Time

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

mean(time.ratio.mean)
sd(time.ratio.mean)

quantile(time.ratio.mean, c(0.025, .975))

mean(time.ratio.mean) - mean(Time.ILEC)/mean(Time.CLEC) # Bias

# Significance
p <- (sum(result >= observed) + 1) / ( N + 1)

ggplot(data.table(result), aes(result)) +
  geom_histogram(aes(fill = ..count..), bins = 30) +
  geom_vline(xintercept = observed, col = "darkorange", lty = 2, lwd = 1.5, alpha = .4) +
  scale_y_continuous(labels = comma) +
  labs(title = "Permutation Test: Difference in means")

