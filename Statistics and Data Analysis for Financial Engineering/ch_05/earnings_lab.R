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

#### Earnings Lab

library("Ecdat")

?CPSch3

data(CPSch3)
dimnames(CPSch3)[[2]]

male.earnings <- CPSch3[CPSch3[, 3] == "male", 2]
sqrt.male.earnings <- sqrt(male.earnings)
log.male.earnings <- log(male.earnings)
one.fourth.male.earnings = male.earnings^(1/4)

par(mfrow = c(2, 2))
qqnorm(male.earnings, datax = T, main = "untransformed")
qqnorm(sqrt.male.earnings, datax = T,
       main = "square-root transformed")
qqnorm(log.male.earnings, datax = T,
       main = "log-transformed")
qqnorm(one.fourth.male.earnings, datax = T,
       main = "1/4 power")

par(mfrow = c(2, 2))
boxplot(male.earnings, datax = T, main = "untransformed")
boxplot(sqrt.male.earnings, datax = T,
        main = "square-root transformed")
boxplot(log.male.earnings, datax = T,
        main = "log-transformed")
boxplot(one.fourth.male.earnings, datax = T,
        main = "1/4 power")

par(mfrow = c(2, 2))
plot(density(male.earnings), main = "untransformed")
plot(density(sqrt.male.earnings),
     main = "square-root transformed")
plot(density(log.male.earnings),
     main = "log-transformed")
plot(density(one.fourth.male.earnings),
     main = "1/4 power")


library(MASS)
par(mfrow = c(1, 1))
boxcox(male.earnings ~ 1)
boxcox(male.earnings ~ 1, lambda = seq(0.3, 0.45, 1 / 100))

bc <- boxcox(male.earnings ~ 1, lambda = seq(0.3, 0.45, by = 1 / 100),
             interp = F)
ind <- (bc$y == max(bc$y))
ind2 <- (bc$y > max(bc$y) - qchisq(0.95, df = 1) / 2)
bc$x[ind]
bc$x[ind2]

plot(bc)

mli = which.max( bc$y ) # what is the largest value of the likelihood
ml_lam = bc$x[mli] # extract the value of lambda that gives the largest likelihood
print(sprintf("Maxlikihood lambda value= %10.6f", ml_lam))
print(range(bc$x[ind2]))

library(car)
plot( density( bcPower( male.earnings, ml_lam ) ) ) # plot the transformed variables


fit_sstd = sstdFit( male.earnings, hessian=T )
