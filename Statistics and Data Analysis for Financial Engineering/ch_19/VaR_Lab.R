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

setwd("D:/Projects/MSDS-RiskAnalytics/Module_07")

library(Ecdat)

data(SP500, package = "Ecdat")

n <- 2783

# non-parametric VaR and ES

SPreturn <- SP500$r500[(n - 999):n]
year <- 1981 + (1:n) * (1991.25 - 1981) / n
year <- year[(n - 999):n]

alpha <- 0.05

q <- as.numeric(quantile(SPreturn, alpha))

VaR_nonp <- -20000 * q

IEVaR <- (SPreturn < q)
sum(IEVaR)

ES_nonp <- -20000 * sum(SPreturn * IEVaR) / sum(IEVaR)

options(digits = 5)
VaR_nonp
ES_nonp

# parametric VaR and ES

library(MASS)

fitt <- fitdistr(SPreturn, "t")
param <- as.numeric(fitt$estimate)

mean <- param[1]
df <- param[3]
sd <- param[2] * sqrt((df) / (df - 2))

lambda <- param[2]
qalpha <- qt(alpha, df = df)
VaR_par <- -20000 * (mean + lambda * qalpha)
es1 <- dt(qalpha, df = df) / (alpha)
es2 <- (df + qalpha^2) / (df - 1)
es3 <- -mean + lambda * es1 * es2
ES_par <- 20000 * es3

VaR_par
ES_par

##################################################
############  Code for Figure 19.2  ##############
##################################################
data(SP500,package="Ecdat")
library("fGarch")
n = 2783
SPreturn = SP500$r500[(n-999):n]
year = 1981 + (1:n)* (1991.25-1981)/n
year = year[(n-999):n]
length(SPreturn)
plot(SPreturn)
par(mfrow=c(2,3))
qqnorm(SPreturn)
n = length(SPreturn)
grid = (1:n)/(n+1)

pdf("SPreturns_tplot.pdf",width=6,height=5)  ###  Figure 19.2
par(mfrow=c(1,1))
qqplot(SPreturn, qt(grid,df=2.9837),main="t-probability plot, 
   df=2.9837",xlab="data",ylab="t-quantiles")
abline(lm(qt(c(.25,.75),df=2.9837)~quantile(SPreturn,c(.25,.75))))
graphics.off()


library(rugarch)
library(fGarch)

garch.t <- ugarchspec(mean.model = list(armaOrder = c(0,0)),
                      variance.model = list(garchOrder = c(1, 1)),
                      distribution.model = "std")

sp.garch.t <- ugarchfit(data = SPreturn, spec = garch.t)

show(sp.garch.t)

pred <- ugarchforecast(sp.garch.t, data = SPreturn, n.ahead = 1)

pred

alpha <- 0.05
nu <- as.numeric(coef(sp.garch.t)[5])
q <- qstd(alpha, mean = fitted(pred), sd = sigma(pred), nu = nu)

VaR <- -20000 * q; VaR

lambda <- sigma(pred) / sqrt( nu / (nu-2) )
qalpha <- qt(alpha, df = nu)
es1 <- dt(qalpha, df = nu) / alpha
es2 <- ( nu + qalpha^2) / ( nu - 1 )
es3 <- -mean + lambda*es1*es2

ES_par <- 20000 * es3; ES_par


####################################################################
############  Code for Example 19.4 and Figure 19.3  ###############
####################################################################
data(SP500,package="Ecdat")
library("fGarch")
n = 2783
alpha = 0.05
SPreturn = SP500$r500[(n-999):n]
year = 1981 + (1:n)* (1991.25-1981)/n
year = year[(n-999):n]
B = 5000
VaRs=matrix(0,nrow=B,ncol=4)
set.seed(38751)
ptm1 <- proc.time()
for (i in (1:B))
{
  returns_b = sample(SPreturn,1000,replace=TRUE)
  q_b = as.numeric(quantile(returns_b,.05))
  VaR_nonp_b = -20000*q_b
  IEVaR_b = (returns_b < q_b)
  ES_nonp_b = -20000 * sum(returns_b*IEVaR_b) / sum(IEVaR_b)
  
  fitt_b = fitdistr(returns_b,"t")
  param_b = as.numeric(fitt_b$estimate)
  mean_b = param_b[1]
  df_b = param_b[3]
  sd_b = param_b[2]*sqrt( (df_b)/(df_b-2) )
  lambda_b = param_b[2]
  qalpha_b = qt(.05,df=df_b)
  VaR_par_b = -20000*(mean_b + lambda_b*qalpha_b)
  es1_b = dt(qalpha_b,df=df_b)/(alpha)
  es2_b = (df_b + qalpha_b^2) / (df_b - 1)
  es3_b = -mean_b+lambda_b*es1_b*es2_b
  ES_par_b = 20000*es3_b
  VaRs[i,]=c(VaR_nonp_b,VaR_par_b,ES_nonp_b,ES_par_b)
}
ptm2 = proc.time()
(ptm2  -  ptm1)/60
colnames(VaRs) = c("VaR_nonp","VaR_par","ES_nonp","ES_par")
write.table(VaRs, file = "SP500VaRBoot.txt",row.names=F) 

quantile(VaRs[,1],c(.025,.975))
quantile(VaRs[,2],c(.025,.975))
quantile(VaRs[,3],c(.025,.975))
quantile(VaRs[,4],c(.025,.975))

####################################################
############  Code for Example 19.5  ###############
####################################################

data(SP500,package="Ecdat")
library("fGarch")
# daily observations from 1981?01 to 1991?04 
# number of observations : 2783 
# daily return S&P500 (change in log index) 
n = 2783
SPreturn = SP500$r500[(n-999):n]
length(SPreturn)
year = 1981 + (1:n)* (1991.25-1981)/n
year = year[(n-999):n]
#  parametric estimation
fitt = fitdistr(SPreturn,"t")
param = as.numeric(fitt$estimate)
mean = param[1]
df = param[3]
sd = param[2]*sqrt( (df)/(df-2) )
alpha = 0.05
n = length(SPreturn)
#  GARCH estimation
library(rugarch)
garch.t = ugarchspec(mean.model=list(armaOrder=c(0,0)),
                     variance.model=list(garchOrder=c(1,1)),
                     distribution.model = "std")
sp.garch.t = ugarchfit(data=SPreturn, spec=garch.t)
show(sp.garch.t)
# plot(sp.garch.t, which = 2)
#pred = as.numeric(predict(fit_garch,n.ahead=1))
pred = ugarchforecast(sp.garch.t, data=SPreturn, n.ahead=1) ; pred
nu = as.numeric(coef(sp.garch.t)[5])
q = qstd(alpha, mean = fitted(pred), 
         sd = sigma(pred), nu = nu )
VaR = -20000*q ; VaR
lambda = sigma(pred)/sqrt( (nu)/(nu-2) )
qalpha = qt(alpha, df=nu)
es1 = dt(qalpha, df=nu)/(alpha)
es2 = (nu + qalpha^2) / (nu - 1)
es3 = -mean+lambda*es1*es2
ES_par = 20000*es3 ; ES_par


pdf("SP500_GARCH_CondSD.pdf",width=6,height=5)   ##  Figure 19.3
par(mfrow=c(1,1), mar=c(5.1, 4.6, 4.1, 2.1))
plot(year,sigma(sp.garch.t),type="l",ylab=expression(hat(sigma)[t]))
points(year[1000]+1/365,sigma(pred),cex=4,pch="*")
abline(h=sd,lty=5,col=2)
legend("topright",c("marginal SD","conditional SD","next day's conditional SD"),
       lty=c(5,1,NA),pch=c(NA,NA,"*"),pt.cex=2.5, col=c(2,1,1))
graphics.off()


######################################################################################
############  Code for Example 19.7 and Figures 19.4, 19.5, and 19.6   ###############
######################################################################################
data(SP500,package="Ecdat")
library("fGarch")
# daily observations from 1981?01 to 1991?04 
# number of observations : 2783 
# daily return S&P500 (change in log index) 
n0 = 2783
SPreturn = SP500$r500[(n0-999):n0]
year = 1981 + (1:n0)* (1991.25-1981)/n0
year = year[(n0-999):n0]
n = length(SPreturn)
x=sort(SPreturn)
plot(SPreturn,xlim=c(0,50),type="b")
pdf("SP500_tail_plots.pdf",width=6,height=5)   ##  Figure 19.4  ##
par(mfrow=c(2,2))
m = 50
xx = log((1:m)/n)
yy = log(-x[1:m])
plot(xx ,yy,  xlab="log(k/n)",
     ylab=expression(paste("log(-",R[(k)],")" )),main="(a) m=50")
fit = lm(yy~xx)
text(-4,-2,paste("slope = ",round(fit$coef[2],3)))
text(-4,-2.25,paste("a = ",round(-1/fit$coef[2],3)))
abline(fit,lwd=2)
m = 100
xx = log((1:m)/n)
yy = log(-x[1:m])
plot(xx ,yy,  xlab="log(k/n)",
     ylab=expression(paste("log(-",R[(k)],")" )),main="(b) m=100")
fit = lm(yy~xx)
text(-3.5,-2,paste("slope = ",round(fit$coef[2],3)))
text(-3.5,-2.35,paste("a = ",round(-1/fit$coef[2],3)))
abline(fit,lwd=2)
m = 200
xx = log((1:m)/n)
yy = log(-x[1:m])
plot(xx ,yy,  xlab="log(k/n)",
     ylab=expression(paste("log(-",R[(k)],")" )),main="(c) m=200")
fit = lm(yy~xx)
text(-2.8,-2,paste("slope = ",round(fit$coef[2],3)))
text(-2.8,-2.35,paste("a = ",round(-1/fit$coef[2],3)))
abline(fit,lwd=2)
m = 300
xx = log((1:m)/n)
yy = log(-x[1:m])
plot(xx ,yy,  xlab="log(k/n)",
     ylab=expression(paste("log(-",R[(k)],")" )),main="(d) m=300")
fit = lm(yy~xx)
text(-3,-2,paste("slope = ",round(fit$coef[2],3)))
text(-3,-2.45,paste("a = ",round(-1/fit$coef[2],3)))
abline(fit,lwd=2)
graphics.off()

q = quantile(x,.1)
print(q,digits=3)
print(-20000*q,digits=3)
a= 1.975
print(1/a,digits=4)
fitt = fitdistr(SPreturn,"t")
param = as.numeric(fitt$estimate)
mean = param[1]
df = param[3]
sd = param[2]*sqrt( (df)/(df-2) )
lambda = param[2]
alpha = seq(.002,.2,.0001)
qalpha = qt(alpha,df=df)
VaR_par = -20000*(mean + lambda*qalpha)
VaR_norm = -20000*(mean(x)+ sd(x)* qnorm(alpha))
pdf("var_sp500_est.pdf",width=6,height=5)                 ##  Figure 19.5  ##
par(mfrow=c(1,1))
plot(alpha, -20000*q * (.1/alpha)^(1/a),type="l",lwd=2,
     xlab=expression(alpha),ylab=expression(paste("VaR(",alpha,")")),
     ylim=c(0,1700),log="x" )
a.hill=2.2
lines(alpha, -20000*q * (.1/alpha)^(1/a.hill),lty=1,lwd=2,col="red")
lines(alpha,VaR_par,lty=1,lwd=2,col="blue")
lines(alpha,VaR_norm,lty=1,lwd=2,col="purple")
legend("topright",c("polynomial tail: regression","polynomial tail: Hill","t","normal"),
       col=c("black", "red", "blue", "purple"),lwd=2,lty=1)
graphics.off()

data(SP500,package="Ecdat")
library("fGarch")
# daily observations from 1981?01 to 1991?04 
# number of observations : 2783 
# daily return S&P500 (change in log index) 
n0 = 2783
SPreturn = SP500$r500[(n0-999):n0]
year = 1981 + (1:n0)* (1991.25-1981)/n0
year = year[(n0-999):n0]
n = length(SPreturn)
x=sort(SPreturn)
plot(SPreturn,type="b")
c = seq(quantile(SPreturn,.025),quantile(SPreturn,.25),by=.001)
nc = 0*c
den=nc
hill=nc
for (i in 1:length(nc))
{
  ind =  (x<c[i]) 
  nc[i] = sum(ind)
  den[i] = sum(  log(x[ind]/c[i])  )
  hill[i] = nc[i]/den[i]
}
pdf("var_sp500_hill.pdf",width=5.25,height=2.25)      ##  Figure 19.6  ##
par(mfrow=c(1,3))
plot(nc,hill,type="b",main="(a)",ylab="Hill estimator")
plot(nc,hill,xlim=c(25,120),ylim=c(2,2.4),type="b",main="(b)",ylab=
       "Hill estimator")
plot(nc,hill,xlim=c(60,100),ylim=c(2,2.4),type="b",main="(c)",ylab=
       "Hill estimator")
hill[8:11]
graphics.off()

qnorm()
