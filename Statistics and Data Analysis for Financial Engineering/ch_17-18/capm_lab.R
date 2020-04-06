##########################################################################################
###########  R script for Chapter 17   ###################################################
###########  of Statistics and Data Analysis for Financial Engineering, 2nd Edition ######
###########  by Ruppert and Matteson  ####################################################
##########################################################################################

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

###  Figures 17.1 and 17.2 were produced by Matlab for the book "Statistics and Finance"  ###

##################################################
############  Code for Figure 17.3   ###############
##################################################

library(Ecdat)
library(quadprog)
data(CRSPday)
R = 100*CRSPday[,4:6]  #  convert to percentages
mean_vect = apply(R,2,mean)
cov_mat = cov(R)
sd_vect = sqrt(diag(cov_mat))
Amat = cbind(rep(1,3),mean_vect)  # set the constraints matrix
muP = seq(.05,.14,length=300)  # set of 300 possible target values
# for the expect portfolio return
sdP = muP # set up storage for std dev's of portfolio returns
weights = matrix(0,nrow=300,ncol=3) # storage for portfolio weights
for (i in 1:length(muP))  # find the optimal portfolios for each target expected return
{
  bvec = c(1,muP[i])  # constraint vector
  result =
    solve.QP(Dmat=2*cov_mat,dvec=rep(0,3),Amat=Amat,bvec=bvec,meq=2)
  sdP[i] = sqrt(result$value)
  weights[i,] = result$solution
}
pdf("sml_derNew.pdf",width=6,height=5)            
plot(sdP,muP,type="l",xlim=c(1.0,1.325),ylim=c(0.076,.12),lty=3, lwd = .5)  #  plot
# the efficient frontier (and inefficient portfolios
# below the min var portfolio)
mufree = 1.3/253 # input value of risk-free interest rate
points(0,mufree,cex=4,pch="*")  # show risk-free asset
sharpe =( muP-mufree)/sdP # compute Sharpe's ratios
ind = (sharpe == max(sharpe)) # Find maximum Sharpe's ratio
options(digits=3)
weights[ind,] #  print the weights of the tangency portfolio
lines(c(0,2),mufree+c(0,2)*(muP[ind]-mufree)/sdP[ind],lwd=4,lty=1, col = "blue")
# show line of optimal portfolios
points(sdP[ind],muP[ind],cex=4,pch="*") # show tangency portfolio
ind2 = (sdP == min(sdP)) # find the minimum variance portfolio
ind3 = (muP > muP[ind2])
lines(sdP[ind3],muP[ind3],type="l",xlim=c(0,.25),
      ylim=c(0,.3),lwd=3, col = "red")  #  plot the efficient frontier
text(sd_vect[1],mean_vect[1],"GE",cex=1.15)
text(sd_vect[2],mean_vect[2],"IBM",cex=1.15)
text(sd_vect[3],mean_vect[3],"Mobil",cex=1.15)
w = seq(0, 1.5,len = 500)
covjm = cov(R[,3], R %*% weights[ind])
mup2 = w * weights[ind] %*% mean_vect + (1 - w) * mean_vect[3]
sdp2 = sqrt(w^2 * sdP[ind]^2 + (1 - w)^2 * sd_vect[3]^2 + 2 * w * (1 - w) * covjm)
lines(sdp2, mup2, lwd = 3, col = "purple")
legend("topleft", c("CML", "efficient frontier","portfolios of tangency and Mobil","tangency portfolio"),
       col = c("blue","red","purple","black"),
       lty = c(1,1,1,NA), lwd = c(4,3,3,NA), pch = c(NA,NA,NA,"*"), pt.cex = c(1,1,1,4)  )
graphics.off()

#####################################################
############  Code for Example 17.3   ###############
#####################################################
dat = read.csv(paste0(data.dir, "capm2.csv"), header=T)
attach(dat)
n = dim(dat)[1]
EX_R_sp500 = Close.sp500[2:n]/Close.sp500[1:(n-1)] - 1  - Close.tbill[2:n]/(100*253) 
EX_R_msft = Close.msft[2:n]/Close.msft[1:(n-1)] - 1  - Close.tbill[2:n]/(100*253) 
fit = lm(EX_R_msft~EX_R_sp500)
options(digits=3)
summary(fit)

fit_NoInt = lm(EX_R_msft~EX_R_sp500-1)
options(digits=3)
summary(fit_NoInt)

##################################################
############  Code for R lab   ###############
##################################################
dat = read.csv(paste0(data.dir, "Stock_Bond_2004_to_2006.csv"),header=T)
prices = dat[,c(5,7,9,11,13,15,17,24)]
n = dim(prices)[1]

dat2 =  as.matrix(cbind(dat[(2:n),3]/365,
                        100*(prices[2:n,]/prices[1:(n-1),] - 1)))
names(dat2)[1] = "treasury"
risk_free = dat2[,1]
ExRet = dat2[,2:9] - risk_free
market = ExRet[,8]
stockExRet = ExRet[,1:7]

fit_reg = lm(stockExRet~market)
summary(fit_reg)
res = residuals(fit_reg)
pairs(res)
options(digits=3)
betas=fit_reg$coeff[2,]




##################################################
############  Section 17.9.1  ####################
##################################################

dat = read.csv(paste0(data.dir, "AlphaBeta.csv"))
alpha = dat$alpha
beta = dat$beta
library(linprog)
B1 = 0.25
B2 = B1
M = length(alpha)
AmatLP1 = diag(1,nrow= 2*M)
AmatLP2 = c(rep(1,M),rep(-1,M))
AmatLP3 = c(beta,-beta)
AmatLP = rbind(AmatLP1,AmatLP2,AmatLP3)
bvecLP = c(rep(B1,M),rep(B2,M),1,0)
cLP =  c(alpha,-alpha)
const.dir = c(rep("<=",2*M),"=","=")
resultLP = solveLP(cvec=cLP, bvec = bvecLP, Amat = AmatLP,
                   lpSolve = TRUE, const.dir = const.dir, maximum=TRUE, verbose = 4)
solution = resultLP$solution
w = solution[1:M] - solution[-(1:M)]
sum(w)
w %*% beta
AmatLP %*% solution
options(digits = 2)
w[1:10]
w[11:20]
w[21:30]
w[31:40]
w[41:50]
w %*% alpha
