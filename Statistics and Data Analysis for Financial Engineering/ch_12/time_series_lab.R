library(Ecdat)
library(tseries)
library(forecast)

setwd("D:/Projects/MSDS-RiskAnalytics/Module_05")

data(Mishkin, package = "Ecdat")

y <- as.vector(Mishkin[, 1])

par(mfrow = c(1,2))
acf(y)
acf(diff(y))

Box.test(diff(y), lag = 10, type = "Ljung-Box")

data(bmw, package = "evir")

Box.test(bmw, lag = 5, type = "Ljung-Box")

fitAR1 <- arima(bmw, order = c(1, 0, 0))
summary(fitAR1)

acf(bmw)

Box.test(residuals(fitAR1), lag = 5, type = "Ljung-Box", fitdf = 1)

fit <- arima(y, order = c(1, 0, 0))
Box.test(fit$residuals, type = "Ljung", lag = 24, fitdf = 1)

library(forecast)

auto.arima(diff(y), max.p = 20, max.q = 0, d = 0, ic = "aic")

auto.arima(diff(y), max.p = 20, max.q = 0, d = 0, ic = "bic")

par(mfrow = c(1,1))
plot.ts(y)

fitMA3 <- arima(diff(y), order = c(0, 0, 3))
fitMA3

set.seed(4631)
y1 <- arima.sim(n = 500, list(ar = c(0.4)))
y2 <- cumsum(y1)
y3 <- cumsum(y2)

par(mfrow = c(3,1))
plot.ts(y1)
plot.ts(y2)
plot.ts(y3)

adf.test(y)
pp.test(y)

kpss.test(y)

?adf.test

x <- rnorm(1000)  # no unit-root
adf.test(x)

y <- diffinv(x)   # contains a unit-root
adf.test(y)

TT <- 100
wn <- rnorm(TT) # white noise
tseries::adf.test(wn)

tseries::adf.test(wn, k=0)

intercept <- 1
wnt <- wn + 1:TT + intercept
tseries::adf.test(wnt)

rw <- cumsum(rnorm(TT))
tseries::adf.test(rw)

tseries::adf.test(rw, k=0)

pdf("inflation_AR1_acf.pdf",width=9,height=5)
#
par(mfrow=c(1,2))
acf(y,main="Inflation rate")
acf(fit$resid,main="Residuals from AR(1)")
#
graphics.off()


x1 = as.vector(ARMAacf(ar=c(0.5,-0.3), lag.max=10))
x2 = as.vector(ARMAacf(ar=c(0.5,0.15), lag.max=10))
x3 = as.vector(ARMAacf(ar=c(0.15,0.8), lag.max=10))

pdf("ar2acf.pdf",width=8,height=6)
#
par(mfrow=c(1,1))
plot(0:10,x1,xlab="lag",ylab="ACF", main= "ACF of three AR(2) processes",cex.axis=1.5,
     cex.lab=1.5,cex=2,cex.main=1.5,pch="*",type="b",ylim=c(-.5,1))
lines(0:10,x2,cex.axis=1.5, cex.lab=1.5,cex=2,pch="o",type="b")
lines(0:10,x3,cex.axis=1.5, cex.lab=1.5,cex=2,pch="x",type="b")

abline(h=0)
legend(6,-.1,c("(0.5, -0.3)", "(0.5, 0.15)","(0.15, 0.8)"), pch=c("*","o","x"), cex=1.5, box.lty=0)
#
graphics.off()

data(Mishkin, package = "Ecdat")
y = as.vector(Mishkin[,1]) 
x = diff(y)
logn = log(length(x))

#####  Fitting AR(p) models  #####  
resultsdiff = matrix(0,nrow=20,ncol=3)
for (i in 1:20){
  fit = arima(x,order=c(i,0,0))
  resultsdiff[i,1] = i  
  resultsdiff[i,2] = fit$aic
  resultsdiff[i,3] = resultsdiff[i,2] + (logn-2)*i
}

pdf("inflation_diff_arfits.pdf",width=8,height=6)
#
plot(resultsdiff[,1],resultsdiff[,2],xlab="p",ylab="criterion",cex.lab=1.35,cex.axis=1.35,
     main="AIC and BIC for AR fits to changes in inflation rate",
     cex.main=1.35,cex=2,pch="*",ylim=c(2440,2560),type='b')
points(resultsdiff[,1],resultsdiff[,3],pch="o",cex=2,type='b')
legend(12,2565,c("AIC","BIC"),pch=c("*","o"),cex=2,box.lty=0)
#
graphics.off()


options(digits=5)

auto.arima(y, max.p = 20, max.q = 0,  d = 0, ic = "aic")
auto.arima(y, max.p = 20, max.q = 0,  d = 0, ic = "bic")



options(digits=7)

auto.arima(diff(y), max.p = 20, max.q = 0, d = 0, ic = "aic")

auto.arima(diff(y), max.p = 20, max.q = 0, d = 0, ic = "bic")


data(Mishkin, package = "Ecdat")
y = as.vector(Mishkin[,1]) 

pdf("inflationAR7_res_acf.pdf", width=7, height=6)
#
par(mfrow=c(2,2))
plot(ts(y, start=1950, frequency=12),ylab="Inflation Rate",type="l",xlab="Year",cex.lab=1.5,
     cex.axis=1.5,cex.main=1.3,main="(a)")
acf(y,main="(b)")
fitAR7 = arima(y,c(7,0,0))
plot(ts(fitAR7$resid, start=c(1950,2), frequency=12),ylab="Inflation Rate",type="l",xlab="Year",cex.lab=1.5,
     cex.axis=1.5,cex.main=1.3,main="(c)")
acf(fit$resid,main="(d)")
#
graphics.off()

data(Mishkin, package = "Ecdat")
y = as.vector(Mishkin[,1]) 
x = diff(y)
logn = log(length(x))

#####  Fitting MA(q) models  #####  
resultsdiff = matrix(0,nrow=9,ncol=3)
for (i in 1:9){
  fit = arima(x,order=c(0,0,i))
  resultsdiff[i,1] = i  
  resultsdiff[i,2] = fit$aic
  resultsdiff[i,3] = resultsdiff[i,2] + (logn-2)*i
}

pdf("inflation_diff_mafits.pdf",width=8,height=6)
#
par(mfrow=c(1,1))
plot(resultsdiff[,1],resultsdiff[,2],xlab="q",ylab="criterion",cex.lab=1.35,cex.axis=1.35,
     main="AIC and BIC for MA fits to changes in inflation rate",
     cex.main=1.35,cex=2,pch="*",ylim=c(2445,2500),type='b')
points(resultsdiff[,1],resultsdiff[,3],pch="o",cex=2,type='b')
legend(1.2,2498,c("AIC","BIC"),pch=c("*","o"),cex=2,box.lty=0)
#
graphics.off()

data(Capm,package="Ecdat")  
rf=Capm$rf
diffrf=diff(rf)
acf(diffrf)
arima(rf,c(2,1,0))

res = matrix(0,nrow=9,ncol=4)
i = 1
for (p in 0:2)
{
  for (q in 0:2)
  {
    res[i,1] = p
    res[i,2] = q
    fit = arima(diffrf,c(p,0,q))
    res[i,4] = AIC(fit,k=logn) +1290
    res[i,3] = AIC(fit) +1290
    i=i+1
  }
}
options(digits=3)
res 

data(Mishkin,package="Ecdat")
y = as.vector(Mishkin[,1]) 

year = seq(1950 + 1/12,1990+11/12,1/12)
n=length(year)
logn=log(n)

fit=arima(y,c(0,1,3))

pred.infl = predict(fit, n.ahead = 100, se.fit = TRUE)
t1 = 300:491
t2 = 492:(492+49+50)
year = seq(1950 + 1/12,2001+61/12,1/12)


pdf("Inflation_predict.pdf",width=7,height=6)
#
plot(year[t1],y[t1],ylim=c(-10,18),type="b",xlim=c(1975,1999),
     xlab="year",ylab="Inflation rate",cex.axis=1.15,cex.lab=1.15)
points(year[t2], pred.infl$pred,type="p",pch="*")
lines(year[t2], pred.infl$pred - 2*pred.infl$se)
lines(year[t2], pred.infl$pred + 2*pred.infl$se)
legend(1975,-3,c("data","predictions","lower CL","upper CL"),cex=1.2,
       box.lty=0,pch=c("o","*",NA,NA),lty=c(NA,NA,1,1))
#
graphics.off()

fit_diff=arima(diff(y),c(0,0,3))

pred.infl_diff =predict(fit_diff, n.ahead = 100, newxreg = NULL,
                        se.fit = TRUE)
t1 = 300:491
t2 = 492:(492+49+50)

pdf("Diff_Inflation_predict.pdf",width=7,height=6)
#
plot(year[t1],diff(y)[t1],xlim=c(1975,1999),ylim=c(-9,15),type="b",
     xlab="year",ylab="Change in inflation rate",cex.axis=1.5,cex.lab=1.5)
points(year[t2], pred.infl_diff$pred,type="p",pch="*")
lines(year[t2], pred.infl_diff$pred - 2*pred.infl_diff$se)
lines(year[t2], pred.infl_diff$pred + 2*pred.infl_diff$se)
legend(1975,14,c("data","predictions","lower CL","upper CL"),cex=1.2,
       box.lty=0,pch=c("o","*",NA,NA),lty=c(NA,NA,1,1))
#
graphics.off()

data(Mishkin,package="Ecdat")
infl = as.vector(Mishkin[,1]) 

year = seq(1950 + 1/12,1990+11/12,1/12)
n=length(year)
logn=log(n)

fit_diff=arima(diff(infl),c(0,0,3))

pred.infl_diff =predict(fit_diff, n.ahead = 100, newxreg = NULL, se.fit = TRUE)
t1 = 300:491
t2 = 492:(492+49+50)

resid = fit_diff$resid[488:490]
coeff = as.vector(fit_diff$coef[1:3])
mu = as.vector(fit_diff$coef[4])
niter = 50000
n.ahead = 30
futureobs = matrix(0,nrow=niter,ncol=n.ahead)
future_int = futureobs

set.seed(1234576)
for (i in 1:niter)
{
  errors = sample(fit_diff$resid, n.ahead, replace = TRUE)
  errors = c(resid,errors)
  for (j in 1:n.ahead)
  {
    futureobs[i,j] = mu + errors[j+3] + errors[j+2]*coeff[1]+ errors[j+1]*coeff[2] + errors[j]*coeff[3]
    if (j > 1)
    {
      future_int[i,j] = future_int[i,j-1] + futureobs[i,j]
    }
    if (j==1){future_int[i,j] = futureobs[i,j]
    }
  }
}
future_mean = apply(futureobs,2,mean)
ul = 0*(1:n.ahead)
ll =ul
for (k in 1:n.ahead)
{
  ul[k] = quantile(futureobs[,k],.975)
  ll[k] = quantile(futureobs[,k],.025)
}

pdf("inflation_forecasts_sim.pdf")
#
plot(1:n.ahead,ul,ylim=c(-10,10),type="b",lwd=2,xlab="month ahead",ylab="rate",cex.axis=1.5,cex.lab=1.5)
lines(ll,type="b",lwd=2)
lines(1:n.ahead, pred.infl_diff$pred[1:n.ahead] - 1.96*pred.infl_diff$se[1:n.ahead],type="b",lty=3)
lines(1:n.ahead, pred.infl_diff$pred[1:n.ahead] + 1.96*pred.infl_diff$se[1:n.ahead],type="b",lty=3)
lines(1:n.ahead, future_mean,lwd=2,lty=2)
#
graphics.off()

pdf("inflation_forecasts_sim_random.pdf")
#
plot(1:n.ahead,futureobs[1,],ylim=c(-12,12),
     type="b",xlab="month ahead",ylab="rate",cex.axis=1.5,cex.lab=1.5,lwd=2)
for (i in 2:5)
{
  lines(1:n.ahead,futureobs[i,],type="b",lty=i,lwd=2)
}
#
graphics.off()

plot(1:n.ahead,future_int[1,],ylim=c(-8,18),
     type="b",xlab="month ahead",ylab="rate",cex.axis=1.5,cex.lab=1.5,lwd=2)
for (i in 2:5)
{
  lines(1:n.ahead,future_int[i,],type="b",lty=i,lwd=2)
}

ul_int = 0*(1:n.ahead)
ll_int =ul_int
for (k in 1:n.ahead)
{
  ul_int[k] = quantile(future_int[,k],.975)
  ll_int[k] = quantile(future_int[,k],.025)
}
future_mean_int = apply(future_int,2,mean)

pdf("inflation_forecasts_sim_integrated.pdf")
#
plot(1:n.ahead,ul_int,ylim=c(-5,15),type="b",lwd=2,xlab="month ahead",ylab="rate",cex.axis=1.5,cex.lab=1.5)
lines(ll_int,type="b",lwd=2)
lines(future_mean_int)
#
graphics.off()


data(bmw,package="evir")

pdf("BMW_pacf.pdf",width=6,height=5)
#
pacf(bmw,main="Sample PACF for daily BMW stock log returns")
#
graphics.off()

data(Mishkin,package="Ecdat")
infl = as.vector(Mishkin[,1])  

pdf("Inflation_pacf.pdf",height=5,width=6)
#
pacf(diff(infl), main = "Change in inflation rate")
#
graphics.off()

