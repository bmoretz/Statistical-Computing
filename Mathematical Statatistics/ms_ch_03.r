library(data.table)
library(ggplot2)

### Beer/Wings Consumption

Beerwings <- data.table(read.csv(paste0(data.dir, "Beerwings.csv"),
                               header = T))

with(Beerwings, {
  tapply(Hotwings, Gender, mean)
})

observed <- 14.5333 - 9.3333 # store observed mean difference

hotwings <- Beerwings$Hotwings

N <- 10^5 - 1

result <- numeric(N)

for( i in 1:N)
{ #sample of size 15, from 1 to 30, without replacement
  index <- sample(30, size = 15, replace = F)
  result[i] <- mean(hotwings[index]) - mean(hotwings[-index])
}

ggplot(data.table(result), aes(result)) +
  geom_histogram(aes(fill = ..count..)) +
  geom_vline(xintercept = observed, col = "darkorange", linetype = 3, lwd = 1.3) +
  scale_y_continuous(labels = comma) +
  labs(title = "Beer/Wings Resampling")

p <- (sum(result >= observed) + 1) / (N + 1) # p-value
v <- p*(1 - p) / N # variance of the samples

### Verizon Repair Times

Verizon <- data.table(read.csv(paste0(data.dir, "Verizon.csv"),
                                 header = T))

with(Verizon, {
  tapply(Time, Group, mean)
})

Time.CLEC <- Verizon[Group == "CLEC"]$Time
Time.ILEC <- Verizon[Group == "ILEC"]$Time
Time <- c(Time.CLEC, Time.ILEC)

observed <- mean(Time.ILEC) - mean(Time.CLEC) # store the observed mean difference

N <- 10e4 - 1

result <- numeric(N)

sample.size <- length(Time.ILEC) + length(Time.CLEC)

for( i in 1:N)
{
  index <- sample(sample.size, length(Time.ILEC), replace = F)
  result[i] <- mean(Time[index]) - mean(Time[-index])
}

ggplot(data.table(result), aes(result)) +
  geom_histogram(aes(fill = ..count..)) +
  geom_vline(xintercept = observed, col = "darkorange", linetype = 3, lwd = 1.3) +
  scale_y_continuous(labels = comma) +
  labs(title = "Verizon Repair Times - Resampling")

p <- ( sum(result <= observed) + 1 ) / ( N + 1 )
v <- p*(1 - p) / N
