library(xts)
library(ggplot2)
library(ggthemes)

theme_set(theme_light())

# Theme Overrides
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             legend.position = "top", legend.title = element_blank())

data.dir <- "D:/Projects/MSDS-RiskAnalytics/datasets/"

setwd("D:/Projects/MSDS-RiskAnalytics/Module_06")

data(Mishkin,package="Ecdat")
infl= as.vector(Mishkin[,1])  # pai1 = one-month inflation rate (in percent, annual rate) 

data(SP500,package="Ecdat")
SPreturn = SP500$r500

library(evir)  # for emplot
data(Garch,package="Ecdat")
attach(Garch)

data(Capm,package="Ecdat")
difflogrf=diff(log(Capm$rf))

diffdm = diff(dm)  #  Deutsch mark
diffbp = diff(bp)  #  British pound
diffcd = diff(cd)  #  Canadian dollar
diffdy = diff(dy)  #  Japanese yen

n1  = length(SPreturn)
year_SP = 1981 + (1:n1)* (1991.25-1981)/n1
d1 = seq(as.Date("1981-01-01"), as.Date("1991-04-01"), by=1)
date_SP = d1[seq(1,length(d1),by=length(d1)/n1)]

n2 = length(diffdm)
year_dm = 1980 + (1:n2)* (1987.5-1980)/n2
d2 = seq(as.Date("1980-01-01"), as.Date("1987-07-01"), by=1)
date_dm = d2[seq(1,length(d2),by=length(d2)/n2)]

n3 = length(difflogrf)
year_rf = 1960 + (1:n3) * (2003 - 1960)/n3
d3 = seq(as.Date("1960-01-01"), as.Date("2003-01-01"), by=1)
date_rf = d3[seq(1,length(d3),by=length(d3)/n3)]

n4 = length(infl)
year_infl = 1950+1/12 + (1:n4) * (1991-1950-1/12)/n4
d4 = seq(as.Date("1950-02-01"), as.Date("1991-01-01"), by=1)
date_infl = d4[seq(1,length(d4),by=length(d4)/n4)]


pdf("garch_examples.pdf", width = 9, height = 7)
#
lwdfact=6
par(mfrow=c(2,2))
plot(xts(abs(SPreturn),date_SP),main="(a) S&P 500 daily return",xlab="year",type="l", 
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,ylab="|log return|",ylim=c(0,.08), 
     minor.ticks = F, major.ticks="year", major.format="%Y")
mod = loess( abs(SPreturn)~year_SP,span=.25)
lines(xts(predict(mod),date_SP),lwd=lwdfact,col=2 )
#
plot(xts(abs(diffbp),date_dm),xlab="year",ylab="|change in rate|",main=
       "(b) BP/dollar exchange rate",type="l",cex.axis=1.5,cex.lab=1.5,cex.main=1.5, 
     minor.ticks = F, major.ticks="year", major.format="%Y")
mod = loess( abs(diffbp)~year_dm,span=.25)
lines(xts(predict(mod),date_dm),lwd=lwdfact,col=2 )
#
plot(xts(abs(difflogrf),date_rf),main="(c) Risk-free interest rate",xlab="year",
     ylab="|change in log(rate)|",type="l",cex.axis=1.5,cex.main=1.5,cex.lab=1.5, 
     minor.ticks = F, major.ticks="year", major.format="%Y")
mod = loess( abs(difflogrf)~year_rf,span=.25)
lines(xts(predict(mod),date_rf),lwd=lwdfact,col=2 )
#
plot(xts(abs(infl-mean(infl)),date_infl),ylab="|rate - mean(rate)|",xlab="year",
     type="l",cex.axis=1.5,cex.main=1.5,cex.lab=1.5,main="(d) Inflation rate", 
     minor.ticks = F, major.ticks="year", major.format="%Y")
mod = loess( abs(infl-mean(infl))~year_infl,span=.3)
lines(xts(predict(mod),date_infl),lwd=lwdfact,col=2 )
#
graphics.off()


##################################################
############  Code for Example 14.1  #############
############  Code for Figure 14.2  ##############
##################################################

n = 10200
e = rnorm(n)
a = e
y = e
sig2 = e^2
omega = 1
alpha = 0.55
phi = 0.8
mu = 0.1
omega/(1-alpha) ; sqrt(omega/(1-alpha))

set.seed("7484")
for (t in 2:n)
{
  a[t] = sqrt(sig2[t])*e[t]
  y[t] = mu + phi*(y[t-1]-mu) + a[t]
  sig2[t+1] = omega + alpha * a[t]^2
}

pdf("garch01.pdf", width = 9, height = 6)
#
par(mfrow=c(2,4))
plot(e[10001:n],type="l",xlab="t",ylab=expression(epsilon),main="(a) white noise")
plot(sqrt(sig2[10001:n]),type="l",xlab="t",ylab=expression(sigma[t]),
     main="(b) conditional std dev")
plot(a[10001:n],type="l",xlab="t",ylab="a",main="(c) ARCH")
plot(y[10001:n],type="l",xlab="t",ylab="y",main="(d) AR+ARCH")
acf(a[10001:n],main="(e) ARCH")
acf(a[10001:n]^2,main="(f) ARCH squared")
acf(y[10001:n],main="(g) AR+ARCH")
acf(y[10001:n]^2,main="(h) AR+ARCH squared")
#
graphics.off()


##################################################
############  Code for Figure 14.3  ##############
##################################################

n = 10500
set.seed("2340")
e = rnorm(n)
a = e
y = e
sig2= e^2
omega = 1
alpha = 0.08
beta = 0.90
phi = 0.8
mu = 0.1

for (t in 2:n)
{
  a[t] = sqrt(sig2[t])*e[t]
  y[t] = mu + phi*(y[t-1]-mu) + a[t]
  sig2[t+1] = omega + alpha * a[t]^2 + beta*sig2[t]
}

pdf("garch03.pdf", width = 9, height = 6)
#
par(mfrow=c(2,4))
plot(e[10001:n],type="l",xlab="t",ylab=expression(epsilon),main="(a) white noise")
plot(sqrt(sig2[10001:n]),type="l",xlab="t",ylab=expression(sigma[t]),
     main="(b) conditional std dev")
plot(a[10001:n],type="l",xlab="t",ylab="a",main="(c) GARCH")
plot(y[10001:n],type="l",xlab="t",ylab="y",main="(d) AR+GARCH")
acf(a[10001:n],main="(e) GARCH")
acf(a[10001:n]^2,main="(f) GARCH squared")
acf(y[10001:n],main="(g) AR+GARCH")
acf(y[10001:n]^2,main="(h) AR+GARCH squared")
#
graphics.off()


##################################################
############  Code for Figure 14.3  ##############
##################################################

n = 10500
set.seed("2340")
e = rnorm(n)
a = e
y = e
sig2= e^2
omega = 1
alpha = 0.08
beta = 0.90
phi = 0.8
mu = 0.1

for (t in 2:n)
{
  a[t] = sqrt(sig2[t])*e[t]
  y[t] = mu + phi*(y[t-1]-mu) + a[t]
  sig2[t+1] = omega + alpha * a[t]^2 + beta*sig2[t]
}

pdf("garch03.pdf", width = 9, height = 6)
#
par(mfrow=c(2,4))
plot(e[10001:n],type="l",xlab="t",ylab=expression(epsilon),main="(a) white noise")
plot(sqrt(sig2[10001:n]),type="l",xlab="t",ylab=expression(sigma[t]),
     main="(b) conditional std dev")
plot(a[10001:n],type="l",xlab="t",ylab="a",main="(c) GARCH")
plot(y[10001:n],type="l",xlab="t",ylab="y",main="(d) AR+GARCH")
acf(a[10001:n],main="(e) GARCH")
acf(a[10001:n]^2,main="(f) GARCH squared")
acf(y[10001:n],main="(g) AR+GARCH")
acf(y[10001:n]^2,main="(h) AR+GARCH squared")
#
graphics.off()



##################################################
############  Code for Example 14.2  #############
############  Code for Figure 14.4  ##############
############  Code for Figure 14.5  ##############
##################################################
library(xts)
library(rugarch)
data(bmw, package="evir")

arma.garch.norm = ugarchspec(mean.model=list(armaOrder=c(1,0)),
                             variance.model=list(garchOrder=c(1,1)))
bmw.garch.norm = ugarchfit(data=bmw, spec=arma.garch.norm)
show(bmw.garch.norm)

length(bmw)

plot(bmw.garch.norm, which="all")

pdf("garch04.pdf", width = 6, height = 9)
#
par(mfrow = c(3,2))
for(i in c(1,3,10,11,8,9)) plot(bmw.garch.norm, which=i)
#
graphics.off()


library(MASS)
e = residuals(bmw.garch.norm, standardize=TRUE)
fitdistr(e,"t")


n = length(e)
grid = (1:n)/(n+1)

pdf("bmwGarch_11_tplot.pdf",width=4,height=4)
#
par(mfrow=c(1,1))
qqplot(sort(as.numeric(e)), qt(grid,df=4),
       main="t-plot, df=4",xlab= "Standardized residual quantiles",
       ylab="t-quantiles")
abline(   lm(   qt(c(.25,.75),df=4)~quantile(e,c(.25,.75))   )   )
#
graphics.off()


arma.garch.t = ugarchspec(mean.model=list(armaOrder=c(1,0)),
                          variance.model=list(garchOrder=c(1,1)),
                          distribution.model = "std")
bmw.garch.t = ugarchfit(data=bmw, spec=arma.garch.t)
show(bmw.garch.t)


#
par(mfrow = c(3,2))
for(i in c(1,3,10,11,8,9)) plot(bmw.garch.t, which=i)
#



##################################################
############  Code for Figure 14.6  ##############
##################################################

rho = function(a,b){ c(1,a*(1-a*b-b^2)/(1-2*a*b-b^2) * (a+b)^(0:9)) }

target = .5

a1 =.1
b1 = uniroot(f=function(b){rho(a1,b)[2] - target},interval = c(0,1-a1))$root

a2=.3
b2 = uniroot(f=function(b){rho(a2,b)[2] - target},interval = c(0,1-a2))$root

a3=.5
b3 = uniroot(f=function(b){rho(a3,b)[2] - target},interval = c(0,1-a3))$root

pdf("garch11ACF.pdf",width=6,height=5)
#
plot(0:10,rho(a1,b1),type="b",ylim=c(0,1),lty=1,lwd=2,
     ylab=expression(paste(rho[a^2],"(h)")),
     xlab="h" )
lines(0:10,rho(a2,b2),type="b",lty=2,lwd=2)
lines(0:10,rho(a3,b3),type="b",lty=3,lwd=2)
legend("topright",c(expression(paste(alpha," = 0.10, ",beta," = 0.894")),
                    expression(paste(alpha," = 0.30, ",beta," = 0.604")),
                    expression(paste(alpha," = 0.50, ",beta," = 0.000")) ),
       lty=1:3,lwd=2)
#
graphics.off()




##################################################
############  Code for Figure 14.7  ##############
##################################################

rho(.10,.86)

data(bmw,package="evir")
res = residuals(arima(bmw,order=c(1,0,0)))

pdf("bmwACFsquared.pdf",width=6,height=5)
#
acf(res^2)
#
graphics.off()

##################################################
############  Code for Figure 14.8  ##############
##################################################

x = seq(-3,3,.005)
gamma = c(-.5,-.2,0,.12,.3,.9)

pdf("leverage_functions.pdf",width=6,height=4)
#
par(mfrow=c(2,3))
for (i in 1:length(gamma)){
  gama = toString(gamma[i])
  plot(x,abs(x) - gamma[i] * x,ylab=expression(paste(g[gamma],"(x)")),
       main=paste("gamma = ",gamma[i]),type="l",lwd=2  )  
}
#
graphics.off()


##################################################
############  Code for Example 14.3  #############
##################################################
library(xts)
library(rugarch)
data(bmw,package="evir")
arma.aparch.t = ugarchspec(mean.model=list(armaOrder=c(1,0)),
                           variance.model=list(model="apARCH", 
                                               garchOrder=c(1,1)),
                           distribution.model = "std")
bmw.aparch.t = ugarchfit(data=bmw, spec=arma.aparch.t,solver="hybrid")
show(bmw.aparch.t)

-5.8985 * 6146 - (-5.8983 * 6146)
-5.9073 * 6146 - (-5.9048 * 6146)

#
par(mfrow = c(3,2))
for(i in c(1,3,10,11,8,9)) plot(bmw.aparch.t, which=i)
#




##################################################
############  Code for Example 14.4  #############
############  Code for Figure 14.9  ##############
##################################################

nelsonplosser = read.csv(paste0(data.dir, "nelsonplosser.csv"), header = TRUE)
nrow(nelsonplosser)
new_np = na.omit(nelsonplosser) 
n = nrow(new_np) ; n
attach(new_np)

fit.lm1 = lm(diff(log(sp)) ~ diff(log(ip)) + diff(bnd))
summary(fit.lm1)
resid.lm1 = resid(fit.lm1)

library(forecast)
auto.arima(resid.lm1)

xregressors = cbind(diff(log(ip)), diff(bnd))
fit.arma = arima(diff(log(sp)), xreg=xregressors, order=c(0,0,1))
print(fit.arma)
resid.arma = fit.arma$res

par(mfrow=c(1,2))
acf(resid.arma)
acf(resid.arma^2)

library(tseries)
?garch
fit.arch = garch(resid.arma, order = c(0,1), trace=F)
summary(fit.arch)
resid.arch = fit.arch$residuals
sigma.arch = as.numeric(fit.arch$fitted.values[,1])

par(mfrow=c(1,2))
ts.plot(sigma.arch)
acf(na.omit(resid.arch))

nelploss.arch.std.res = as.numeric(resid.arch/sigma.arch)[-1] 

pdf("nelsonPlosser_acf.pdf",width=6,height=5)
#
par(mfrow=c(2,2))
acf(rstudent(fit.lm1),main="(a) regression: residuals")
acf(rstudent(fit.lm1)^2,main="(b) regression: squared residuals")
acf(nelploss.arch.std.res,main="(c) MA/ARCH: residuals")
acf(nelploss.arch.std.res^2,main="(d) MA/ARCH: squared residuals")
#
graphics.off()

fit.lm3 = lm(diff(log(sp)) ~ diff(log(ip)) + diff(bnd),
             weights = 1/sigma.arch^2) 
summary(fit.lm3)

par(mfrow=c(1,1))
plot(fitted(fit.lm1)[-1], fitted(fit.lm3)) ; abline(0,1)


##################################################
############  Code for Example 14.5  #############
############  Code for Figure 14.10  #############
##################################################

library(xts)
library(rugarch)
data(bmw, package="evir")
n = length(bmw) ; n
date = attr(bmw,"time")
#year = 1973 + (1996 + 7/12 - 1973)*(1:n)/n
origin1 = 4100
date[origin1]
date[origin2]
nahead = 1500

garch.norm = ugarchspec(mean.model=list(armaOrder=c(0,0)),
                        variance.model=list(garchOrder=c(1,1)))

bmw.garch.norm.1 = ugarchfit(data=bmw[1:origin1], spec=garch.norm)
pred1 = ugarchforecast(bmw.garch.norm.1, data = bmw[1:origin1], n.ahead = nahead)
head(fitted(pred1))
head(sigma(pred1))
plot(pred1, which = 1)

bmw.garch.norm.2 = ugarchfit(data=bmw[1:origin2], spec=garch.norm)
pred2 = ugarchforecast(bmw.garch.norm.2, data = bmw[1:origin2], n.ahead = nahead)
head(fitted(pred2))
head(sigma(pred2))
plot(pred2)
library(xts)
BMW = xts(bmw, attr(bmw,"times"))

pdf("bmw_garchForecast.pdf",height=5,width=6)
#
lwd1 = 4
plot(BMW,type="l", xlab="year",ylab="log return", ylim=c(-.13,.21), xlim=c(date[3393],date[5219]), main="Forecasting BMW returns") #,xlim=c(1986,1992),ylim=c(-.13,.21),
#
lines(date[(origin2+1):(origin2+nahead)],fitted(pred2)+1.96*sigma(pred2),col=5,lwd=lwd1)
lines(date[(origin2+1):(origin2+nahead)],fitted(pred2)-1.96*sigma(pred2),col=5,lwd=lwd1)
#
lines(date[(origin1+1):(origin1+nahead)],fitted(pred1)+1.96*sigma(pred1),col=2,lwd=lwd1,lty=2)
lines(date[(origin1+1):(origin1+nahead)],fitted(pred1)-1.96*sigma(pred1),col=2,lwd=lwd1,lty=2)
#
legend("topleft",c("11-15-87","9-18-88"),lty=c(1,2),lwd=lwd1,col=c(5,2))
#
graphics.off()




##################################################
############  Code for Figure 14.11  #############
##################################################

data(CRSPday, package="Ecdat")
CRSPday = ts(CRSPday, start = c(1989, 1), frequency = 253) 
ibm   = CRSPday[,5] * 100
crsp  = CRSPday[,7] * 100
Y = cbind(ibm, crsp) 

pdf("mvData.pdf",width=9,height=8)
#
par(mfrow = c(2,1))
plot(Y[,1], type = 'l', xlab = "year", ylab = "return (%)", main = "(a)")
plot(Y[,2], type = 'l', xlab = "year", ylab = "return (%)", main = "(b)")
#
graphics.off()





##################################################
############  Code for Figure 14.12  #############
##################################################

pdf("acfccf.pdf",width=9,height=8)
#
layout(rbind(c(1,2), c(3,3)),widths=c(1,1,2),heights=c(1,1))
acf(as.numeric(Y[,1]), ylim = c(-0.1,0.1), main = "(a)")
acf(as.numeric(Y[,2]), ylim = c(-0.1,0.1), main = "(b)")
ccf(as.numeric(Y[,1]),as.numeric(Y[,2]), 
    type = c("correlation"), main = "(c)", ylab = "CCF", lag = 20)
#
graphics.off()

cor(ibm, crsp)


# MV Box-Ljung test for bivariate returns
source("SDAFE2.R")
mLjungBox(Y, 5)


# Fitting an VAR(1) model (for simplicity)
#Y = Y - colMeans(Y)   # (mean zero returns)
fit.AR1 = ar(Y, aic = FALSE, order.max = 1)
fit.AR1

A = fit.AR1$resid[-1,]
mLjungBox(A, 5)

#mLjungBox(A, 5, df.adj = 2*2*1)
#acf(A, ylim = c(-0.1,0.2))


##################################################
############  Code for Figure 14.13  #############
##################################################

pdf("acfSq.pdf",width=7,height=7) #
#
par(mfrow = c(2,2))
acf(A[,1]^2, ylim = c(-0.05,0.35), main = "(a)")
acf(A[,2]^2, ylim = c(-0.05,0.35), main = "(b)")
ccf(A[,1]^2,A[,2]^2, type = c("correlation"), main = "(c)", ylab = "CCF", lag = 20)
acf(A[,1]*A[,2], ylim = c(-0.05,0.35), main = "(d)")
#
graphics.off()



##################################################
############  Code for Figure 14.14  #############
##################################################

#########  EWMA model  ############
source("SDAFE2.R")
EWMA.param = est.ewma(lambda.0 = 0.95, innov=A) 
EWMA.param$lambda.hat
EWMA.Sigma = sigma.ewma(lambda = EWMA.param$lambda.hat, innov=A)

pdf("EWMAfit.pdf",width=9,height=8)
#
par(mfrow = c(2,2))
plot(ts(EWMA.Sigma[1,1,]^.5, start = c(1989, 1), frequency = 253), 
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(a) ", {hat(sigma)[{"1,t"}]})))
plot(ts(EWMA.Sigma[1,2,], start = c(1989, 1), frequency = 253),  
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(b) ", {hat(sigma)[{"12,t"}]})))
plot(ts(EWMA.Sigma[1,2,]/(sqrt(EWMA.Sigma[1,1,]* EWMA.Sigma[2,2,])), 
        start = c(1989, 1), frequency = 253),  
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(c) ", {hat(rho)[{"12,t"}]})))
points(ts(mvwindow.cor(A[,1],A[,2], win = 126)$correlation, 
          start = c(1989, 1), frequency = 253), 
       type = 'l', col = 2, lty = 2, lwd=2)
plot(ts(EWMA.Sigma[2,2,]^.5, start = c(1989, 1), frequency = 253),  
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(d) ", {hat(sigma)[{"2,t"}]})))
#
graphics.off()


################################################################
######  MGARCH via DOC assuming GARCH(1,1)-t components ########
######  Includes Orthogonal (PCA) GARCH(1,1) component model ###
################################################################
#library(mnormt)
library(fGarch)
source("SDAFE2.R")
DOC.fit = doc.garch(E = A, L = 4., c = 2.25, theta.ini = NULL)

names(DOC.fit)

DOC.fit$coef.pca
DOC.fit$coef.doc

DOC.fit$W.hat
DOC.fit$U.hat
DOC.fit$M.hat
solve(DOC.fit$U.hat)

acf(DOC.fit$Z.hat^2)
acf(DOC.fit$S.hat^2)


##################################################
############  Code for Figure 14.15  #############
##################################################

pdf("PCAfit.pdf",width=9,height=8)
#
par(mfrow = c(2,2))
plot(ts(DOC.fit$Sigma.pca[1,1,]^.5, start=c(1989,1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(a) ", {hat(sigma)[{"1,t"}]})))
plot(ts(DOC.fit$Sigma.pca[2,1,], start=c(1989,1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(b) ", {hat(sigma)[{"12,t"}]})))
plot(ts(DOC.fit$Sigma.pca[2,1,]/(sqrt(DOC.fit$Sigma.pca[1,1,]*
                                        DOC.fit$Sigma.pca[2,2,])), 
        start=c(1989,1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, ylim = c(0.1,0.9),
     main = expression(paste("(c) ", {hat(rho)[{"12,t"}]})))
points(ts(mvwindow.cor(A[,1],A[,2], win = 126)$correlation, 
          start = c(1989, 1), frequency = 253), 
       type = 'l', col = 2, lty = 2, lwd = 2)
plot(ts(DOC.fit$Sigma.pca[2,2,]^.5, start=c(1989,1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(d) ", {hat(sigma)[{"2,t"}]})))
#
graphics.off()


##################################################
############  Code for Figure 14.16  #############
##################################################

pdf("DOCfit.pdf",width=9,height=8)
#
par(mfrow = c(2,2))
plot(ts(DOC.fit$Sigma.doc[1,1,]^.5, start=c(1989,1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(a) ", {hat(sigma)[{"1,t"}]})))
plot(ts(DOC.fit$Sigma.doc[2,1,], start=c(1989,1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(b) ", {hat(sigma)[{"12,t"}]})))
plot(ts(DOC.fit$Sigma.doc[2,1,]/(sqrt(DOC.fit$Sigma.doc[1,1,]*
                                        DOC.fit$Sigma.doc[2,2,])), 
        start=c(1989,1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, ylim = c(0.1,0.9),
     main = expression(paste("(c) ", {hat(rho)[{"12,t"}]})))
points(ts(mvwindow.cor(A[,1],A[,2], win = 126)$correlation, 
          start=c(1989,1), frequency=253), 
       type = 'l', col = 2, lty = 2,lwd=2)
plot(ts(DOC.fit$Sigma.doc[2,2,]^.5, start=c(1989,1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(d) ", {hat(sigma)[{"2,t"}]})))
#
graphics.off()


acf(DOC.fit$S.hat^2)

# original obs squared
DOC.test(A^2, m = 5)
# After PCA rotation 
DOC.test(DOC.fit$Z.hat^2, m = 5)
# After DOC rotation
DOC.test(DOC.fit$S.hat^2, m = 5)


##################################################
############  Code for Figure 14.17  #############
##################################################

#########  DCC_E model  ###########
library(fGarch)
source("SDAFE2.R")
DCCe.fit = fit.DCCe(theta.0=0.95, innov=A)

DCCe.fit$coef
DCCe.fit$lambda

pdf("DCCefit.pdf",width=9,height=8)
#
par(mfrow = c(2,2))
plot(ts(DCCe.fit$Sigma.t[1,1,]^.5, start=c(1989, 1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(a) ", {hat(sigma)[{"1,t"}]})))
plot(ts(DCCe.fit$Sigma.t[2,1,], start=c(1989, 1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(b) ", {hat(sigma)[{"12,t"}]})))
plot(ts(DCCe.fit$R.t[2,1,], start=c(1989, 1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(c) ", {hat(rho)[{"12,t"}]})))
points(ts(mvwindow.cor(A[,1],A[,2], win = 126)$correlation, 
          start=c(1989, 1), frequency=253), 
       type = 'l', col = 2, lty = 2, lwd=2)
plot(ts(DCCe.fit$Sigma.t[2,2,]^.5, start=c(1989, 1), frequency=253), 
     type = 'l', xlab = "year", ylab = NULL, 
     main = expression(paste("(d) ", {hat(sigma)[{"2,t"}]})))
#
graphics.off()


###################################
###  Residual/Model Checking  #####
###################################
source("SDAFE2.R")
mLjungBox(A^2, lag=5)

n = dim(A)[1]
d = dim(A)[2]

stdResid.EWMA = matrix(0,n,d)
stdResid.PCA  = matrix(0,n,d)
stdResid.DOC  = matrix(0,n,d)
stdResid.DCCe = matrix(0,n,d)

for(t in 1:n){
  stdResid.EWMA[t,] = A[t,] %*% matrix.sqrt.inv(EWMA.Sigma[,,t])
  stdResid.PCA[t,]  = A[t,] %*% matrix.sqrt.inv(DOC.fit$Sigma.pca[,,t])
  stdResid.DOC[t,]  = A[t,] %*% matrix.sqrt.inv(DOC.fit$Sigma.doc[,,t])
  stdResid.DCCe[t,] = A[t,] %*% matrix.sqrt.inv(DCCe.fit$Sigma.t[,,t])
}

mLjungBox(stdResid.EWMA^2, lag=5)
mLjungBox(stdResid.PCA^2, lag=5)
mLjungBox(stdResid.DOC^2, lag=5)
mLjungBox(stdResid.DCCe^2, lag=5)

mLjungBox(stdResid.EWMA[,1] * stdResid.EWMA[,2], lag=5)
mLjungBox(stdResid.PCA[,1] * stdResid.PCA[,2], lag=5)
mLjungBox(stdResid.DOC[,1] * stdResid.DOC[,2], lag=5)
mLjungBox(stdResid.DCCe[,1] * stdResid.DCCe[,2], lag=5)

