library(data.table)
library(ggplot)
library(ggtheme)
library(Ecdat)
library(faraway)
library(fGarch)

theme_set(theme_light())

setwd("D:/Projects/MSDS-RiskAnalytics/Module_03")

# Theme Overrides
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             legend.position = "top", legend.title = element_blank())

data(SP500, package = "Ecdat")

SPreturn <- SP500$r500
n <- length(SPreturn)
year_SP <- 1981 + (1:n) * (1991.25 - 1981) / n

# Visualize S&P500 log returns
plot(year_SP, SPreturn, main = "S&P 500 daily log returns",
     xlab = "year", type = "l", ylab = "log return")

ggplot(data.table(date = year_SP, return = SPreturn, direction = ifelse(SPreturn >= 0, "red", "green")), 
       aes(date, return, color = direction)) +
  geom_bar(stat = 'identity') +
  theme(legend.position = "none")

set.seed("991155")

edf_norm <- ecdf(rnorm(150))
pdf("normalcdfplot.pdf", width = 6, height = 5)
par(mfrow = c(1, 1))
plot(edf_norm, verticals = T, do.p = F, main = "EDF and CDF")
tt = seq(from = -3, to = 3, by = 0.01)
lines(tt, pnorm(tt),lty = 2, lwd = 2, col = "red")
legend(1.5, 0.2, c("EDF", "CDF"), lty = c(1, 2),
       lwd = c(1.5, 2), col = c("black", "red"))
graphics.off()

dm <- Garch$dm
diffdm <- diff(dm) # Ducsch mar
pdf("dm_halfnormal.pdf", width = 7, height = 6)
halfnorm(abs(diffdm), main = "changes in DM/dollar exchange rate",
         ylab = "Sorted data")
graphics.off()


rf <- Capm$rf
diffrf <- diff(rf) # Risk-free rate
pdf("rf_halfnormal.pdf", width = 7, height = 6)
halfnorm(abs(diffrf), main = "changes in Risk-Free rate",
         ylab = "Sorted data")
graphics.off()


qqplot(SPreturn, diffdm, xxlab = "S&P return",
       ylab = "change in DM/dollar rate", main = "(a)")
xx <- quantile(SPreturn, c(0.25, 0.75))
yy <- quantile(diffdm, c(0.25, 0.75))
slope <- (yy[2] - yy[1]) / (xx[2] - xx[1])
inter <- yy[1] - slope*xx[1]
abline(inter, slope, lwd = 2)


shapiro.test(SPreturn)

EuStockMarkets

data(EuStockMarkets)
mode(EuStockMarkets)
class(EuStockMarkets)
plot(EuStockMarkets)

pdf("EuStocks.pdf", width = 6, height = 5)
plot(EuStockMarkets)
graphics.off()


ggplot(EuStockMarkets, aes(DAX, fill = ..count..)) +
  geom_histogram()

logR <- diff(log(EuStockMarkets))
plot(logR)

plot(as.data.frame(logR))

par(mfrow = c(2, 2))
for(i in colnames(logR)) {
  
  qqnorm(logR[, i], datax = T, main = i)
  qqline(logR[, i], datax = T)
  print(shapiro.test(logR[, i]))
}

n <- dim(logR)[1]
q_grid <- (1:n) / (n+1)
df_grid <- c(1, 4, 6, 10, 20, 30)
index.names <- dimnames(logR)[[2]]

for(i in 1:4) {
  # dev.new()
  par(mfrow = c(3,2))
  for(df in df_grid) {
    qqplot(logR[, i], qt(q_grid, df),
           main = paste(index.names[i], "df = ", df) )
    abline(lm(qt(c(0.25, 0.75), df = df) ~
                quantile(logR[, i], c(0.25, 0.75))))
  }
}


x <- seq(-0.1, 0.1, by = 0.001)
par(mfrow = c(1, 1))
df <- 5
mad_t <- mad(logR[, 1],
             constant = sqrt(df / ( df - 2)) / qt(0.75, df))

plot(density(logR[, 1]), lwd = 2, ylim = c(0, .5), xlim = c(-0.1, -0.05))

lines(x, dstd(x, mean = mean(logR[, 1]), sd = mad_t, nu = df),
      lty = 5, lwd =2, col = "red")
lines(x, dnorm(x, mean = mean(logR[, 1]), sd = sd(logR[, 1])),
      lty = 3, lwd = 4, col = "blue")
legend("topleft", c("KDE", paste("t: df = ", df), "normal"),
       lwd = c(2, 2, 4), lty = c(1, 5, 3),
       col = c("black", "red", "blue"))


# P 8: Try to fit a t-distribution to the CAC index returns
#
# First use QQ plots to estimate a good value for the degrees of freedom:
#

n = dim(logR)[1]
q.grid = (1:n)/(n+1)
df = c(1,4,6,10,20,30)
col_index=3 # the CAC index returns
par(mfrow=c(3,2))
for(j in 1:6){
  qqplot(logR[,col_index], qt(q.grid,df=df[j]), main=paste(index.names[col_index], ", df=", df[j]))
  abline(lm( qt(c(0.25,0.75),df=df[j]) ~ quantile(logR[,col_index],c(0.25,0.75)) ))
}
par(mfrow=c(1,1))

# With that value compare the parametric fit with that dof to the true density:
#

x = seq(-0.1,+0.1,by=0.001)
col_indx = 3 # study the CAC index returns
df = 6
mad_t = mad(logR[, 1], constant=sqrt( df / (df - 2) ) / qt(0.75, df) )
par(mfrow=c(1, 3))

plot( density(logR[,col_indx]), lwd=2, ylim=c(0,60), main='both tails' )
lines(x, dstd(x, mean=mean(logR[,col_indx]), sd=mad_t, nu=df),lty=5, lwd=2, col='red')
lines(x, dnorm(x, mean=mean(logR[,col_indx]), sd=sd(logR[,col_indx])), lty=3, lwd=4, col='blue')
legend( 'topleft', c('KDE', paste('t: df= ', df),'normal'), lwd=c(2,2,4), lty=c(1,5,3), col=c('black', 'red', 'blue'))

plot( density(logR[,col_indx]), lwd=2, ylim=c(0,20), xlim=c(-0.05,-0.0), xlab='', main='left tail' )
lines(x, dstd(x, mean=mean(logR[,col_indx]), sd=mad_t, nu=df),lty=5, lwd=2, col='red')
lines(x, dnorm(x, mean=mean(logR[,col_indx]), sd=sd(logR[,col_indx])), lty=3, lwd=4, col='blue')
plot( density(logR[,col_indx]), lwd=2, ylim=c(0,20), xlim=c(0.0, 0.05), xlab='', main='right tail' )
lines(x, dstd(x, mean=mean(logR[,col_indx]), sd=mad_t, nu=df),lty=5, lwd=2, col='red')
lines(x, dnorm(x, mean=mean(logR[,col_indx]), sd=sd(logR[,col_indx])), lty=3, lwd=4, col='blue')
par(mfrow=c(1, 1))


### McDonnald's Prices and Returns

path.data <- "D:/Projects/MSDS-RiskAnalytics/datasets"
setwd(path.data)

# 4.10.2
# Data Analytics

data <- read.csv("MCD_PriceDaily.csv", header = T)
head(data)

nrow(data)
data[1,]$Date
data[1177,]$Date

adjPrice <- data[, 7]
plot(adjPrice, type = "l", lwd = 2)

n <- nrow(adjPrice)

LogRet = diff(log(adjPrice))

getRetVsT(LogRet)

plot(data$Date[-1], LogRet, type='l', lwd=2, xlab='Date', ylab='logreturn')
grid()

par(mfrow=c(1, 2))
hist(LogRet, 80, freq=FALSE)
qqnorm(LogRet, datax=T)
qqline(LogRet, datax=T)
par(mfrow=c(1, 1))