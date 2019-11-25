library(data.table)
library(ggplot)
library(ggtheme)
library(Ecdat)
library(faraway)
library(fGarch)
library(sn)

theme_set(theme_light())

setwd("D:/Projects/MSDS-RiskAnalytics/datasets")

# Theme Overrides
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             legend.position = "top", legend.title = element_blank())

#### Maximum Likelihood Lab

data(Capm, package = "Ecdat")
x <- diff(Capm$rf)
fitdistr(x, "t")

# classical t-distribution

n <- length(x)
start <- c(mean(x), sd(x), 5)
loglik_t <- function(beta) sum( - dt((x - beta[1]) / beta[2],
                                     beta[3], log = T) + log(beta[2]))
fit_t <- optim(start, loglik_t, hessian = T,
               method = "L-BFGS-B", lower = c(-1, 0.001, 1))
AIC_t <- 2 * fit_t$value + 2 * 3
BIC_t <- 2 * fit_t$value + log(n) * 3
sd_t <- sqrt(diag(solve(fit_t$hessian)))

fit_t$par
sd_t
AIC_t
BIC_t

# standardized t-distribution

loglik_std <- function(beta) sum( - dstd(x, mean = beta[1],
                                         sd = beta[2], nu = beta[3], log = T))
fit_std <- optim(start, loglik_std, hessian = T,
                 method = "L-BFGS-B", lower = c(-0.1, 0.01, 2.1))
AIC_std <- 2 * fit_std$value + 2 * 3
BIC_std <- 2 * fit_std$value + log(n) * 3
sd_std <- sqrt(diag(solve(fit_std$hessian)))

fit_std$par
sd_std
AIC_std
BIC_std

# F-S skewed t-distribution

loglik_sstd <- function(beta) sum( - dsstd(x, mean = beta[1],
                                           sd = beta[2], nu = beta[3], xi = beta[4], log = T))
start <- c(mean(x), sd(x), 5, 1)
fit_sstd <- optim(start, loglik_sstd, hessian = T,
                  method = "L-BFGS-B", lower = c(-0.1, 0.01, 2.1, -2))
AIC_sstd <- 2 * fit_sstd$value + 2 * 3
BIC_sstd <- 2 * fit_sstd$value + log(n) * 3
sd_sstd <- sqrt(diag(solve(fit_sstd$hessian)))

fit_sstd$par
sd_sstd
AIC_sstd
BIC_sstd

# Fitting distributions with MLE

dat <- read.csv("FlowData.csv")
dat <- dat/10000

par(mfrow = c(3,2))

fit_mle <- function(data, desc) {
  x <- data
  x1 <- sort(x)
  fit1 <- sn.mple(y = x1, x = as.matrix(rep(1, length(x1))))
  est1 <- cp2dp(fit1$cp, family = "SN")
  
  plot(x1, dsn(x1, dp = est1),
       type = "l", lwd = 2, xlab = "flow",
       ylab = paste(desc," density"))
  d <- density(x1)
  lines(d$x, d$y, lty = 2, lwd = 2)
  legend(40, 0.034, c("t-model", "KDE"), lty = c(1, 2),
         lwd = c(2, 2))
  n <- length(x1)
  u = (1:n) / (n + 1)
  
  plot(x1, qsn(u, dp = est1), xlab = "data",
       ylab = "skew-t quantiles", main = desc)
  lmfit <- lm(qsn(c(0.25, 0.75), dp = est1) ~ quantile(x1, 
                                                       c(0.25, 0.75)))
  abline(lmfit)
}

fit_mle(dat$Flow1, "Flow 1")
fit_mle(dat$Flow2, "Flow 2")
fit_mle(dat$Flow3, "Flow 3")


# Profile Likelihood

library(MASS)

adj = 1.5
par(mfrow = c(3, 3))
x <- dat$Flow1
x1 <- sort(x)

bcfit1 <- boxcox(x1 ~ 1, lambda = seq(2.6, 4.5, 1 / 100),
                 xlab = expression(alpha))
text(3, -1898.75, "Flow 1")
plot(density((x1^3.5 - 1) / 3.5, adjust = adj), main = "Flow 1")
qqnorm((x1^3.5- 1) / 3.5, datax = T, main = "Flow 1")

# KDE / TKDE

y <- diff(Capm$rf)
diffrf <- y
x1 <- pstd(y, mean = 0.001, sd = 0.0725, nu = 3.34)
x <- qnorm(x1)

par(mfrow = c(1, 1))

d1 <- density(diffrf)
plot(d1$x, d1$y, type = "l", xlab = "y", ylab = "Density (y)",
     lwd = 2)

d2 <- density(x)
ginvx <- qstd(pnorm(d2$x), mean = 0.001, sd = 0.0725, nu = 3.34)
gprime_num <- dstd(ginvx, mean = 0.001, sd = 0.0725, nu = 3.34)
gprime_den <- dnorm(qnorm(pstd(ginvx, mean = 0.001,
                               sd = 0.0725, nu = 3.34)))
gprime <- gprime_num / gprime_den
lines(ginvx, d2$y * gprime, type = "l", lty = 2, col = "red", lwd = 2)
legend("topleft", c("KDE", "TKDE"), lty = c(1, 2), lwd = 2,
       col = c("black", "red"))

### Fitting the Multivariate t-Dist by MLE

library(mnormt)
library(MASS)

par(mfrow = c(1,1))

data(CRSPday, package = "Ecdat")

dat <- CRSPday[, 4:7]
df <- seq(5.25, 6.75, 0.01)

n <- length(df)
loglik <- rep(0, n)
for(i in 1:n){
  fit <- cov.trob(dat, nu = df)
  loglik[i] <- sum(log(dmt(dat, mean = fit$center,
                           S = fit$cov, df = df[i])))
}

aic_t <- -max(2 * loglik) + 2 * (4 + 10 + 1) + 64000
z1 <- ( 2 * loglik > 2 * max(loglik) - qchisq(0.95, 1) )

vals <- data.table( df = df, v_hat = 2 * loglik - 64000)

plot(df, 2 * loglik - 64000, type = "l", cex.axis = 1.5,
     cex.lab = 1.5, ylab = "2 * loglikelihood - 64,000", lwd = 2)
abline(h = 2 * max(loglik) - qchisq(0.95, 1) - 64000)
abline(h = 2 * max(loglik) - 64000)
abline(v = (df[16] + df[17]) / 2)
abline(v = (df[130] + df[131]) / 2)

v_hat <- vals[which(vals$v_hat == max(2*loglik - 64000))]

cor(dat)
c$cor

c <-  cov.trob(dat, nu = v_hat$df)
round(c$center, 7)

fit

### Multivariate Skewed t-Distribution

library(sn)

fit <- mst.mple(y = dat, penalty = NULL)
aic_skewt <- -2 * fit$logL + 64000 + 2 * ( 4 + 10 + 4 + 1)
dp2cp(fit$dp, "st")
aic_skewt

yboot <- dat[sample((1:n), n, replace = T),]

colMeans(yboot)
