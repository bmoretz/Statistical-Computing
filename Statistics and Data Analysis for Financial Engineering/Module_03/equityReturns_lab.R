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

#### Equity Returns Lab

berndtInvest <- read.csv("berndtInvest.csv")
Berndt <- as.matrix(berndtInvest[, 2:5])

cov(Berndt)
cor(Berndt)

pairs(Berndt)

library(MASS)
library(mnormt)

df <- seq(2.5, 8, 0.1)
n <- length(df)

loglik_profile <- rep(0, n)
for(i in 1:n)
{
  fit <- cov.trob(Berndt, nu = df[i])
  mu <-  as.vector(fit$center)
  sigma <- matrix(fit$cov, nrow = 4)
  loglik_profile[i] <- sum(log(dmt(Berndt, mean = fit$center,
                                   S = fit$cov, df = df[i])))
}

aic_t <- -max(2 * loglik_profile) + 2 * (4 + 10 + 1)
z1 <- ( 2 * loglik_profile > 2 * max(loglik_profile) - qchisq(0.95, 1) )
plot(df, 2 * loglik_profile, type = "l", cex.axis = 1.5,
     cex.lab = 1.5, ylab = "2 * loglikelihood - 64,000", lwd = 2)
abline(h = 2 * max(loglik_profile) - qchisq(0.95, 1) - 64000)
abline(h = 2 * max(loglik_profile) - 64000)
abline(v = (df[16] + df[17]) / 2)
abline(v = (df[130] + df[131]) / 2)