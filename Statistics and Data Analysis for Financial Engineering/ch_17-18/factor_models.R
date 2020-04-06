library(data.table)
library(xts)
library(ggplot2)
library(ggthemes)

theme_set(theme_light())

theme_set(theme_light())

# Theme Overrides
theme_update(axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "darkgreen"),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             plot.subtitle = element_text(face = "bold", size = 8, colour = "darkred"),
             legend.title = element_text(size = 12, color = "darkred", face = "bold"),
             legend.position = "right", legend.title.align=0.5,
             panel.border = element_rect(linetype = "solid", 
                                         colour = "lightgray"), 
             plot.margin = unit(c( 0.1, 0.1, 0.1, 0.1), "inches"))

data.dir <- "D:/Projects/MSDS-RiskAnalytics/datasets/"

setwd("D:/Projects/MSDS-RiskAnalytics/Module_09")

###########  R script for Chapter 18  ####################################################
###########  of Statistics and Data Analysis for Financial Engineering, 2nd Edition ######
###########  by Ruppert and Matteson  ####################################################

####################################################################
############  Code for Example 18.2 and Figure 18.1  ###############
####################################################################

datNoOmit = read.table(paste0(data.dir, "treasury_yields.txt"),header=T)
diffdatNoOmit = diff(as.matrix(datNoOmit[,2:12]))
dat=na.omit(datNoOmit)
diffdat = na.omit(diffdatNoOmit)
n = dim(diffdat)[1]
options(digits=5)
pca = prcomp(diffdat)
summary(pca)
pdf("yields01.pdf",width=7,height=6)   ###  Figure 18.1
par(mfrow=c(2,2))
time = c(1/12,.25,.5,1, 2, 3, 5, 7, 10, 20, 30)
plot(time,as.vector(dat[1,2:12]),ylim=c(0,6),type="b",lty=1,lwd=2,
     ylab="Yield",xlab="T",main="(a)") #,log="x",xaxs="r")
lines(time,as.vector(dat[486,2:12]),type="b",lty=2,lwd=2,col="red")
lines(time,as.vector(dat[n+2,2:12]),type="b",lty=3,lwd=2,col="blue")
legend("bottomright",c("07/31/01","07/02/07","10/31/08"),lty=c(1,2,3),lwd=2,
       cex=1, col=c("black","red","blue"))
plot(pca,main="(b)")
plot(time,pca$rotation[,1],,ylim=c(-.8,.8),type="b",lwd=2,ylab="PC",xlab="T",
     main="(c)")
lines(time,pca$rotation[,2],lty=2,type="b",lwd=2,col="red")
lines(time,pca$rotation[,3],lty=3,type="b",lwd=2,col="blue")
lines(0:30,0*(0:30),lwd=1)
legend("bottomright",c("PC 1","PC 2","PC 3"),lty=c(1,2,3),lwd=2,col=c("black","red","blue"))
plot(time,pca$rotation[,1],ylim=c(-.8,.8),type="b",lwd=2,ylab="PC",xlab="T",
     xlim=c(0,3),main="(d)")
lines(time,pca$rotation[,2],lty=2,type="b",lwd=2,col="red")
lines(time,pca$rotation[,3],lty=3,type="b",lwd=2,col="blue")
lines(0:30,0*(0:30),lwd=1)
legend("bottomright",c("PC 1","PC 2","PC 3"),lty=c(1,2,3),lwd=2,col=c("black","red","blue"))
graphics.off()

##################################################
############  Code for Figure 18.2  ###############
##################################################
pdf("yields02.pdf",width=7,height=6)
mu = apply(dat[,2:12],2,mean)
par(mfrow=c(2,2))
plot(time,mu,ylim=c(2.75,4.65),type="b",lwd=4,xlab="T",ylab="Yield",
     main="(a)",xlim=c(0,7),col="black")
lines(time,mu+pca$rotation[,1],lty=5,type="b",lwd=2,col="red")
lines(time,mu-pca$rotation[,1],lty=5,type="b",lwd=2,col="blue")
legend("bottomright",c("mean","mean + PC1",
                       "mean -  PC1"),lty=c(1,5,5),lwd=c(4,2,2),col=c("black","red","blue"))
plot(time,mu,ylim=c(2.75,4.65),type="b",lwd=4,xlab="T",ylab="Yield",
     main="(b)", xlim=c(0,7))
lines(time,mu+pca$rotation[,2],lty=5,type="b",lwd=2,col="red")
lines(time,mu-pca$rotation[,2],lty=5,type="b",lwd=2,col="blue")
legend("bottomright",c("mean","mean + PC2","mean -  PC2"),lty=c(1,5,5),
       lwd=c(4,2,2),col=c("black","red","blue"))
plot(time,mu,ylim=c(2.75,4.65),type="b",lwd=4,xlab="T",ylab="Yield",
     main="(c)",xlim=c(0,7))
lines(time,mu+pca$rotation[,3],lty=5,type="b",lwd=2,col="red")
lines(time,mu-pca$rotation[,3],lty=5,type="b",lwd=2,col="blue")
legend("bottomright",c("mean","mean + PC3",
                       "mean -  PC3"),lty=c(1,5,5),lwd=c(4,2,2),col=c("black","red","blue"))
par(lwd=1)
plot(time,pca$rotation[,4],ylim=c(-.7,.7),type="b",lwd=2,ylab="PC",xlab="T",
     xlim=c(0,30),main="(d)")
lines(time,pca$rotation[,5],lty=5,type="b",lwd=2,col="red")
lines(0:30,0*(0:30),lwd=1)
legend("topright",c("PC 4","PC 5"),lty=c(1,5,5),lwd=2,col=c("black","red"))
graphics.off()


###########################################################
############  Code for Figure 18.3 and 18.4 ###############
###########################################################
pdf("yields_pca_tsplots.pdf",width=7,height=3)  ##  Figure 18.3
par(mfrow=c(1,3))
for (i in 1:3){
  plot(pca$x[,i],main=paste("PC",toString(i)),xlab="day",
       ylab="")
}
graphics.off()

pdf("yields_pca_acfplots2.pdf",width=5,height=5)  ##  Figure 18.4
acf(pca$x[,1:3],ylab="",xlab="lag")
graphics.off()

#####################################################################
############  Code for Example 18.3 and Figure 18.5   ###############
#####################################################################S
equityFunds = read.csv(paste0(data.dir,"equityFunds.csv"))
equityFunds[1:10,]
pairs(equityFunds[,2:9])
pcaEq = prcomp(equityFunds[,2:9])
summary(pcaEq)

pdf("equityFunds_scree.pdf",width=7.2,height=3.5)  ###  Figure 18.5
par(mfrow=c(1,2))
plot(pcaEq,main="(a)")
Names = names(equityFunds)[2:9]
plot(pcaEq$rotation[,1],type="b",ylab="PC",lwd=2,ylim=c(-1.4,2),main="(b)")
lines(pcaEq$rotation[,2],type="b",lty=2,lwd=2,col="red")
lines(pcaEq$rotation[,3],type="b",lty=3,lwd=2,col="blue")
lines(0:9,0*(0:9))
legend("top",c("PC1","PC2","PC3"),lty=c(1,2,3),lwd=2,cex=.65,col=c("black", "red", "blue"))
text(4.35,-1.25, "   EASTEU   LATAM   CHINA   INDIA   ENERGY   MINING   GOLD   WATER",cex=.38)
graphics.off()

#####################################################
############  Code for Example 18.4   ###############
#####################################################
DowJones30 = read.csv(paste0(data.dir, "DowJones30.csv"))
pcaDJ = prcomp(DowJones30[,2:31])
summary(pcaDJ)

#####################################################################
############  Code for Example 18.5 and Figure 18.6   ###############
#####################################################################
pdf("macrofactors.pdf",width=6,height=7)
CPI.dat = read.csv(paste0(data.dir, "CPI.dat.csv"))
IP.dat = read.csv(paste0(data.dir, "IP.dat.csv"))
berndtInvest = read.csv(paste0(data.dir, "berndtInvest.csv"))
berndt = as.matrix(berndtInvest[,-1])   #  1978-01-01 to 1987-12-01
CPI.dat = read.csv(paste0(data.dir, "CPI.dat.csv"))
IP.dat = read.csv(paste0(data.dir, "IP.dat.csv"))
berndt = as.matrix(berndtInvest[,-1])   #  1978-01-01 to 1987-12-01
CPI2 = as.matrix(CPI.dat$CPI[775:900]) #  1977-07-30  to 1987-12-31
CPI = as.data.frame(diff(log(CPI2)))  
names(CPI)[1]="CPI"
IP2 = as.matrix(IP.dat$IP)[703:828,]   #  1977-07-28 to 1987-12-28 
IP = as.data.frame(diff(log(IP2)))
names(IP)[1] = "IP"
CPI_IP = cbind(CPI,IP)
arFit = ar(cbind(CPI,IP))
res = arFit$resid[6:125,]
lmfit = lm(berndt[,2:10]~res[,1]+res[,2])
slmfit = summary(lmfit)
rsq=rep(0,9)
for (i in 1:9){rsq[i]= slmfit[[i]][[8]]}
beta_CPI = lmfit$coef[2,]
beta_IP = lmfit$coef[3,]
par(mfrow=c(1,3))
barplot(rsq,horiz=T,names=names(beta_CPI),main="R squared")
barplot(beta_CPI,hori=T,main="beta CPI")
barplot(beta_IP,hori=T,main="beta IP")
graphics.off()

#########################################################################################
############  Code for Examples 18.6 and 18.7 and Figures 18.7 and 18.8   ###############
#########################################################################################S
#  Uses monthly data from Jan-69 to Dec-98

FF_data = read.table(paste0(data.dir, "FamaFrench_mon_69_98.txt"),header=T)
attach(FF_data)
library("Ecdat")
library("robust")
data(CRSPmon)

ge = 100*CRSPmon[,1] - RF
ibm = 100*CRSPmon[,2] - RF
mobil = 100*CRSPmon[,3] - RF
stocks=cbind(ge,ibm,mobil)
fit = lm(cbind(ge,ibm,mobil)~Mkt.RF+SMB+HML)
options(digits=3)
fit

pdf("FamaFrenchPairs.pdf",width=8,height=8)  ##  Figure 18.7
pairs(cbind(ge,ibm,mobil,Mkt.RF,SMB,HML))
graphics.off()


cor(fit$residuals)
covRob(fit$residuals,cor=T)
cor.test(fit$residuals[,1], fit$residuals[,2])
cor.test(fit$residuals[,1], fit$residuals[,3])
cor.test(fit$residuals[,2], fit$residuals[,3])

pdf("FamaFrenchResidualsPairs.pdf",width=6,height=5)  ## Figure 18.8
pairs(fit$residuals)
graphics.off()


################################################################
########## Code for Example 18.7   #############################
################################################################
sigF = as.matrix(var(cbind(Mkt.RF,SMB,HML)))
bbeta = as.matrix(fit$coef)
bbeta = t( bbeta[-1,])
n=dim(CRSPmon)[1]
sigeps = (n-1)/(n-4) * as.matrix((var(as.matrix(fit$resid))))
sigeps = diag(as.matrix(sigeps))
sigeps = diag(sigeps,nrow=3)
cov_equities = bbeta %*% sigF %*% t(bbeta) + sigeps

options(digits=5)
sigF
bbeta
sigeps
bbeta %*% sigF %*% t(bbeta)
cov_equities
cov(stocks)

###############################################################################
############  Code for Example 18.8 and Figures 18.9 and 18.10  ###############
###############################################################################

berndtInvest = read.csv(paste0(data.dir, "berndtInvest.csv"))
returns = berndtInvest[,-c(1,11,18)]
ind_codes = as.factor(c(3,3,2,1,1,2,3,3,1,2,2,3,1,2,3))
codes = as.matrix(model.matrix(~ind_codes))
codes[,1] =  1 - codes[,2] - codes[,3]
betas = as.data.frame(codes)
colnames(betas) = c("tech","oil","other")
rownames(betas) = colnames(berndtInvest[,-c(1,11,18)])
betas

factors = matrix(0,nrow=120,ncol=3)
for (i in 1:120)
{
  return_data = cbind(t(returns[i,]),betas)
  colnames(return_data)[1] = "return"
  lmfit = lm(return~betas[,1] + betas[,2], data=return_data)
  factors[i,]= lmfit$coef
}

pdf("berndt_cross_section_factors.pdf",width=7,height=3) ###  Figure 18.9
par(mfrow=c(1,3),cex.axis=1.08,cex.lab=1.08,cex.main=1.05)
plot(factors[,1],type="b",lty="dotted",
     lwd=2,xlab="month",ylab="factor",main="market")
plot(factors[,2],lty="dotted",lwd=2,type="b",
     xlab="month",ylab="factor",main="technology")
plot(factors[,3],lty="dotted",lwd=2,type="b",
     xlab="month",ylab="factor",main="oil")
graphics.off()

colnames(factors) = c("market","tech","oil")
cor(factors)
options(digits=2)
sqrt(diag(cov(factors)))

pdf("berndt_cross_section_factors_acf.pdf",width=6,height=6)  ###  Figure 18.10
acf(factors,ylab="",xlab="lag")
graphics.off()

################################################################
############  Code for Examples 18.9 and 18.10   ###############
################################################################

equityFunds = read.csv(paste0(data.dir, "equityFunds.csv"))
fa_none = factanal(equityFunds[,2:9],4,rotation="none")
fa_vari = factanal(equityFunds[,2:9],4,rotation="varimax")
print(fa_none,cutoff=0.1)
print(fa_vari,cutoff=0.1,sort=T)
print(fa_vari,cutoff=0.1)
sum(fa_vari$loadings[,1]^2)
B=fa_vari$loadings[,]

B_none = fa_none$loadings[,]
BB_none = B_none %*% t(B_none)
D_none = diag(fa_none$unique)
Sigma_R_hat = BB_none + D_none
t(B_none) %*% solve(D_none) %*% B_none

options(digits=3)

diff = Sigma_R_hat - fa_vari$corr
product = solve(Sigma_R_hat) %*% fa_vari$corr
max(diff)
min(diff)
eig_diff = eigen(diff)
sort(eig_diff$values)

w = matrix(1/8,nrow =1,ncol=8)
w %*% BB_none %*% t(w)
w %*% fa_vari$corr %*% t(w)
w %*% eig_diff$vectors %*% t(w)

t(B) %*% diag(1/fa_vari$unique) %*% B
t(BB_none) %*% diag(1/fa_none$unique) %*% BB_none

##################################################
############  Code for the R lab   ###############
##################################################

################################################################
########## Code for Section 18.8.1   ###########################
################################################################

yieldDat = read.table("yields.txt",header=T)
maturity = c((0:5),5.5,6.5,7.5,8.5,9.5)
pairs(yieldDat)
par(mfrow=c(4,3))
for (i in 0:11)
{
  plot(maturity,yieldDat[100*i+1,],type="b")
}

eig = eigen(cov(yieldDat))
eig$values
eig$vectors
par(mfrow=c(1,1))
barplot(eig$values)

par(mfrow=c(2,2))
plot(eig$vector[,1],ylim=c(-.7,.7),type="b")
abline(h=0)
plot(eig$vector[,2],ylim=c(-.7,.7),type="b")
abline(h=0)
plot(eig$vector[,3],ylim=c(-.7,.7),type="b")
abline(h=0)
plot(eig$vector[,4],ylim=c(-.7,.7),type="b")
abline(h=0)


library("tseries")
adf.test(yieldDat[,1])

n=dim(yieldDat)[1]
delta_yield = yieldDat[-1,] - yieldDat[-n,]

par(mfrow=c(1,1))
pca_del = princomp(delta_yield)
names(pca_del)
summary(pca_del)
plot(pca_del)

fa_none_5 = factanal(equityFunds[,2:9],5,rotation="none")  ## Does not run


################################################################
########## Code for Section 18.8.2    ##########################
################################################################

#  Extracts daily data 2004-2005 from FamaFrenchDaily.txt

stocks = read.csv("Stock_FX_Bond_2004_to_2005.csv",header=T)
attach(stocks)
stocks_subset=as.data.frame(cbind(GM_AC,F_AC,UTX_AC,MRK_AC))
FF_data0 = read.table("FamaFrenchDaily.txt",header=T)
dat = FF_data0[FF_data0$date > 20040000,]
FF_data = dat[dat$date<20060000,]
FF_data = FF_data[-1,] # delete first row since stocks_diff
# lost a row due to differencing

stocks_diff = as.data.frame(100*apply(log(stocks_subset),2,diff) - FF_data$RF)
names(stocks_diff) = c("GM","Ford","UTX","Merck")
detach(stocks)

fit1 = lm(as.matrix(stocks_diff)~FF_data$Mkt.RF)
summary(fit1)

fit2 = lm(as.matrix(stocks_diff)~FF_data$Mkt.RF +
            FF_data$SMB + FF_data$HML)
summary(fit2)


################################################################
########## Code for Section 18.8.3   ###########################
################################################################

dat = read.csv("Stock_FX_Bond.csv")
stocks_ac = dat[,c(3,5,7,9,11,13,15,17)]
n = length(stocks_ac[,1])
stocks_returns = log(stocks_ac[-1,] / stocks_ac[-n,])
fact = factanal(stocks_returns,factors=2,,rotation="none")
print(fact)
print(fact,cutoff = 0.01)
print(fact,cutoff = 0)

loadings = matrix(as.numeric(loadings(fact)),ncol=2)
unique = as.numeric(fact$unique)
loadings
unique
