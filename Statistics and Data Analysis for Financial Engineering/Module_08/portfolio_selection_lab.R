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

setwd("D:/Projects/MSDS-RiskAnalytics/Module_08")

################################################################
########## Code for figure 16.1   ##############################
################################################################

mu1 = 0.14
mu2 = 0.08
sig1 = 0.2
sig2 = 0.15
rho = 0
rf = 0.06
w = seq(0, 1, len = 500)
means = 0.08 + 0.06 * w
var = sig1^2 * w^2 + sig2^2 * (1 - w)^2
risk = sqrt(var)
ind = !(risk > min(risk))
ind2 = (means > means[ind])
wt = 0.693
meant = 0.08 + 0.06 * wt
riskt = sqrt(sig1^2 * wt^2 + sig2^2 * (1 - wt)^2)

wp = 0.475
meanp = 0.08 + 0.06 * wp
riskp = sqrt(sig1^2 * wp^2 + sig2^2 * (1 - wp)^2)


pdf("portfolioNew.pdf", width = 6, height = 5)

plot(risk, means, type = "l", lwd = 1, xlim = c(0, 0.21), ylim = c(0.0575, 0.145))
abline(v = 0)
lines(risk[ind2], means[ind2], type = "l", lwd = 5, col = "red")
lines( c(0, riskt), c(rf, meant), col = "blue", lwd =2)
lines(c(0,riskp), c(rf,meanp), col = "purple", lwd = 2, lty = 2)
text(riskt, meant, "T", cex = 1.2)
text(sig1, mu1, "R1", cex = 1.2)
text(sig2, mu2, "R2", cex = 1.2)
text(0, rf, "F", cex = 1.2)
text(riskp, meanp, "P", cex = 1.2)
text(min(risk), means[ind], "MV", cex = 1.2)

graphics.off()


################################################################
########## Code for figure 16.2  ###############################
################################################################
pdf("portfolio01New.pdf", width = 6, height = 5)

par(mfrow = c(2, 2))
mu1 = 0.14
mu2 = 0.08
sig1 = 0.2
sig2 = 0.15
rf = 0.06
w = seq(-2, 3, len = 500)
means = 0.08 + 0.06 * w

rho = 0.5
V1 = mu1 - rf
V2 = mu2 - rf
z1 = V1 * sig2^2 - V2 * rho * sig1 *sig2
z2 = V2 * sig1^2 - V1 * rho * sig1 * sig2
wt = z1/(z1 + z2)
meant = 0.08 + 0.06 * wt
riskt = sqrt(sig1^2 * wt^2 + sig2^2 * (1 - wt)^2 + 2 * rho * sig1 * sig2 * wt * (1 - wt))
risk  = sqrt(sig1^2 * w^2 + sig2^2 *   (1 - w)^2 + 2 * rho * sig1 * sig2 * w * (1 - w))
ind = !(risk > min(risk))
ind2 = (means > means[ind])
plot(risk, means, type = "l", lwd = 1, xlim = c(0, 0.3), ylim = c(0.0575, 0.18),
     main = paste("rho = ", rho))
abline(v = 0)
lines(risk[ind2], means[ind2], type = "l", lwd = 5, col = "red")
lines( c(0, riskt), c(rf, meant), col = "blue", lwd =2)
text(riskt, meant, "T", cex = 1.2)

rho = 0.3
V1 = mu1 - rf
V2 = mu2 - rf
z1 = V1 * sig2^2 - V2 * rho * sig1 *sig2
z2 = V2 * sig1^2 - V1 * rho * sig1 * sig2
wt = z1/(z1 + z2)
meant = 0.08 + 0.06 * wt
riskt = sqrt(sig1^2 * wt^2 + sig2^2 * (1 - wt)^2 + 2 * rho * sig1 * sig2 * wt * (1 - wt))
risk  = sqrt(sig1^2 * w^2 + sig2^2 *   (1 - w)^2 + 2 * rho * sig1 * sig2 * w * (1 - w))
ind = !(risk > min(risk))
ind2 = (means > means[ind])
plot(risk, means, type = "l", lwd = 1, xlim = c(0, 0.3), ylim = c(0.0575, 0.18),
     main = paste("rho = ", rho))
abline(v = 0)
lines(risk[ind2], means[ind2], type = "l", lwd = 5, col = "red")
lines( c(0, riskt), c(rf, meant), col = "blue", lwd =2)
text(riskt, meant, "T", cex = 1.2)

rho = 0.0
V1 = mu1 - rf
V2 = mu2 - rf
z1 = V1 * sig2^2 - V2 * rho * sig1 *sig2
z2 = V2 * sig1^2 - V1 * rho * sig1 * sig2
wt = z1/(z1 + z2)
meant = 0.08 + 0.06 * wt
riskt = sqrt(sig1^2 * wt^2 + sig2^2 * (1 - wt)^2 + 2 * rho * sig1 * sig2 * wt * (1 - wt))
risk  = sqrt(sig1^2 * w^2 + sig2^2 *   (1 - w)^2 + 2 * rho * sig1 * sig2 * w * (1 - w))
ind = !(risk > min(risk))
ind2 = (means > means[ind])
plot(risk, means, type = "l", lwd = 1, xlim = c(0, 0.3), ylim = c(0.0575, 0.18),
     main = paste("rho = ", rho))
abline(v = 0)
lines(risk[ind2], means[ind2], type = "l", lwd = 5, col = "red")
lines( c(0, riskt), c(rf, meant), col = "blue", lwd =2)
text(riskt, meant, "T", cex = 1.2)

rho = -0.7
V1 = mu1 - rf
V2 = mu2 - rf
z1 = V1 * sig2^2 - V2 * rho * sig1 *sig2
z2 = V2 * sig1^2 - V1 * rho * sig1 * sig2
wt = z1/(z1 + z2)
meant = 0.08 + 0.06 * wt
riskt = sqrt(sig1^2 * wt^2 + sig2^2 * (1 - wt)^2 + 2 * rho * sig1 * sig2 * wt * (1 - wt))
risk  = sqrt(sig1^2 * w^2 + sig2^2 *   (1 - w)^2 + 2 * rho * sig1 * sig2 * w * (1 - w))
ind = !(risk > min(risk))
ind2 = (means > means[ind])
plot(risk, means, type = "l", lwd = 1, xlim = c(0, 0.3), ylim = c(0.0575, 0.18),
     main = paste("rho = ", rho))
abline(v = 0)
lines(risk[ind2], means[ind2], type = "l", lwd = 5, col = "red")
lines( c(0, riskt), c(rf, meant), col = "blue", lwd =2)
text(riskt, meant, "T", cex = 1.2)
graphics.off()




################################################################
########## Code for Example 16.6 and Figure 16.3  ##############
################################################################
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

bvec = c(1,muP[1])

pdf("quad_prog_plot.pdf",width=6,height=5)              
######## Figure 16.3  #########
par(mfrow = c(1,1))
plot(sdP,muP,type="l",xlim=c(0,2.5),ylim=c(0,.15),lty=3)  #  plot
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
points(sdP[ind2],muP[ind2],cex=2,pch="+") # show min var portfolio
ind3 = (muP > muP[ind2])
lines(sdP[ind3],muP[ind3],type="l",xlim=c(0,.25),
      ylim=c(0,.3),lwd=3, col = "red")  #  plot the efficient frontier
text(sd_vect[1],mean_vect[1],"GE",cex=1.15)
text(sd_vect[2],mean_vect[2],"IBM",cex=1.15)
text(sd_vect[3],mean_vect[3],"Mobil",cex=1.15)

graphics.off()


################################################################
########## Code for Example 16.7 and Figure 16.4  ##############
################################################################
library(Ecdat)
library(quadprog)
data(CRSPday)
R = 100*CRSPday[,4:6]  #  convert to percentages
mean_vect = apply(R,2,mean)
cov_mat = cov(R)
sd_vect = sqrt(diag(cov_mat))
Amat = cbind(rep(1,3),mean_vect,diag(1,nrow=3))  # set the constraints matrix
muP = seq(min(mean_vect)+.0001,max(mean_vect)-.0001,length=300)  
# set of 300 possible target values 
# for the expect portfolio return
sdP = muP # set up storage for standard deviations of portfolio returns
weights = matrix(0,nrow=300,ncol=3) # storage for portfolio weights

for (i in 1:length(muP))  # find the optimal portfolios for each target expected return
{
  bvec = c(1,muP[i],rep(0,3))
  result = 
    solve.QP(Dmat=2*cov_mat,dvec=rep(0,3),Amat=Amat,bvec=bvec,meq=2)
  sdP[i] = sqrt(result$value)
  weights[i,] = result$solution
}

pdf("quad_prog_plotNoShort.pdf",width=6,height=5)       ##########  Figure 16.4  ############

par(mfrow = c(1,1))
plot(sdP,muP,type="l",xlim=c(0,2.5),ylim=c(0,.15),lty=3)  #  plot 
# the efficient frontier (and inefficient frontier)
mufree = 1.3/253 # input value of risk-free interest rate
points(0,mufree,cex=4,pch="*")  # show risk-free asset
sharpe =( muP-mufree)/sdP # compute Sharpe ratios
ind = (sharpe == max(sharpe)) # Find maximum Sharpe ratio
weights[ind,] # Find tangency portfolio
lines(c(0,sdP[ind]),c(mufree,muP[ind]),lwd=4,lty=1, col = "blue") # show line of optimal portfolios
points(sdP[ind],muP[ind],cex=4,pch="*") # show tangency portfolio
ind2 = (sdP == min(sdP)) # find the minimum variance portfolio
points(sdP[ind2],muP[ind2],cex=2,pch="+") # show minimum variance portfolio
ind3 = (muP > muP[ind2])
lines(sdP[ind3],muP[ind3],type="l",xlim=c(0,.25),
      ylim=c(0,.3),lwd=2, col = "red")  #  plot the efficient frontier
text(sd_vect[1],mean_vect[1],"GE",cex=1.15)
text(sd_vect[2],mean_vect[2],"IBM",cex=1.15)
text(sd_vect[3],mean_vect[3],"Mobil",cex=1.51)

graphics.off()


################################################################
###### Code for Examples 16.8 and 16.9 and Figure 16.5  ########
################################################################
library(quadprog)

dat = read.table(paste0(data.dir, "countries.txt"),header =T)
dat = dat[,4:13]

n = dim(dat)[1]
N = dim(dat)[2]
R = 100*(dat[2:n,]/dat[1:(n-1),] - 1)

###########  Short Sales Allowed  #########
mufree = 1/24
mean_vect_TRUE = apply(R,2,mean)
cov_mat_TRUE = cov(R)
nboot = 250
out = matrix(1,nrow=nboot,ncol=2)
mean_out = matrix(1,nrow = nboot,ncol = dim(dat)[2])
set.seed(998877)
for (iboot in (1:nboot))
{
  un = ceiling((n-1)*runif(n-1))
  Rboot = R[un,]
  mean_vect = apply(Rboot,2,mean)
  mean_out[iboot,] = mean_vect
  cov_mat = cov(Rboot)
  sd_vect = sqrt(diag(cov_mat))
  Amat = cbind(rep(1,N),mean_vect) 
  muP = seq(0,2.5,length=300)                              
  sdP = muP 
  weights = matrix(0,nrow=300,ncol=N) 
  for (i in 1:length(muP))  
  {
    bvec = c(1,muP[i])  
    result = 
      solve.QP(Dmat=2*cov_mat,dvec=rep(0,N),Amat=Amat,bvec=bvec,meq=2)
    sdP[i] = sqrt(result$value)
    weights[i,] = result$solution
  } 
  sharpe =( muP-mufree)/sdP 
  ind = (sharpe == max(sharpe)) 
  out[iboot,1] = sharpe[ind]
  wT = weights[ind,]
  sharpe_TRUE = (wT %*% mean_vect_TRUE - mufree) /
    sqrt(wT %*% cov_mat_TRUE %*% wT)
  out[iboot,2] = sharpe_TRUE
}
out_Short = out
gp = cbind(rep("estimated",nboot),rep("actual",nboot))
country_names = c(
  "Hong Kong",
  "Singapore",
  "Brazil",
  "Argentina",
  "UK",
  "Germany",
  "Canada",
  "France",
  "Japan",
  "US")
pdf("countries_boot.pdf",width=7,height=4)      ########  Figure 11.5  ########  
par(mfrow=c(1,2))
boxplot(out_Short~gp,main="(a) Short Sales Allowed",ylim=c(0,.7))
abline(h=.3681,lwd=3,lty=2)

###################  No short sales  ##################

set.seed(998877)
for (iboot in (1:nboot))
{
  un = ceiling((n-1)*runif(n-1))
  Rboot = R[un,]
  mean_vect = apply(Rboot,2,mean)
  cov_mat = cov(Rboot)
  sd_vect = sqrt(diag(cov_mat))
  Amat = cbind(rep(1,N),mean_vect,diag(1,N)) 
  muP = seq(min(mean_vect)+.001,max(mean_vect)-.001,length=300)                              
  sdP = muP 
  weights = matrix(0,nrow=300,ncol=N) 
  for (i in 1:length(muP))  
  {
    bvec = c(1,muP[i],rep(0,N))  
    result = 
      solve.QP(Dmat=2*cov_mat,dvec=rep(0,N),Amat=Amat,bvec=bvec,meq=2)
    sdP[i] = sqrt(result$value)
    weights[i,] = result$solution
  } 
  sharpe =( muP-mufree)/sdP 
  ind = (sharpe == max(sharpe)) 
  out[iboot,1] = sharpe[ind]
  wT = weights[ind,]
  sharpe_TRUE = (wT %*% mean_vect_TRUE - mufree) /
    sqrt(wT %*% cov_mat_TRUE %*% wT)
  out[iboot,2] = sharpe_TRUE
}
out_NoShort = out
boxplot(out_NoShort~gp,main="(b) No Short Sales",ylim=c(0,.7))  ###  Figure 16.5 continued
abline(h=.3503,lwd=3,lty=2)
graphics.off()


################################################################
########## Code for Table 16.1   ###############################
################################################################
###  Because of Monte Carlo variability, the results in the table will not be
###  reproduced exactly.  Increasing nboot will reduce Monte Carlo variability
nboot = 50000
library(bootstrap)
for (i in 1:10){
  bsamp = bootstrap(R[,i],nboot,mean)$thetastar
  print(as.numeric(quantile(bsamp,c(0.025,0.975))))
}


#####################  Example 16.10 and Figure 11.6  #############################

library(quadprog)
dat = read.table("countries.txt",header =T)
dat = dat[,4:13]
n = dim(dat)[1]
N = dim(dat)[2]
R = 100*(dat[2:n,]/dat[1:(n-1),] - 1)
mufree = 1/24
mean_vect_TRUE = apply(R,2,mean)
cov_mat_TRUE = cov(R)
nboot = 250
gp = cbind(rep("estimated",nboot),rep("actual",nboot))
country_names = c(
  "Hong Kong",
  "Singapore",
  "Brazil",
  "Argentina",
  "UK",
  "Germany",
  "Canada",
  "France",
  "Japan",
  "US")

########################  No short sales and no shrinkage  ####################

out = matrix(1,nrow=nboot,ncol=2)
set.seed(998877)
for (iboot in (1:nboot))
{
  un = ceiling((n-1)*runif(n-1))
  Rboot = R[un,]
  mean_vect = apply(Rboot,2,mean)
  cov_mat = cov(Rboot)
  sd_vect = sqrt(diag(cov_mat))
  Amat = cbind(rep(1,N),mean_vect,diag(1,N)) 
  muP = seq(min(mean_vect)+.001,max(mean_vect)-.001,length=300)                              
  sdP = muP 
  weights = matrix(0,nrow=300,ncol=N) 
  for (i in 1:length(muP))  
  {
    bvec = c(1,muP[i],rep(0,N))  
    result = 
      solve.QP(Dmat=2*cov_mat,dvec=rep(0,N),Amat=Amat,bvec=bvec,meq=2)
    sdP[i] = sqrt(result$value)
    weights[i,] = result$solution
  } 
  sharpe =( muP-mufree)/sdP 
  ind = (sharpe == max(sharpe)) 
  out[iboot,1] = sharpe[ind]
  wT = weights[ind,]
  sharpe_TRUE = (wT %*% mean_vect_TRUE - mufree) /
    sqrt(wT %*% cov_mat_TRUE %*% wT)
  out[iboot,2] = sharpe_TRUE
}
out_NoShort = out

##########################  No short sales but with shrinkage  ########################

out = matrix(1,nrow=nboot,ncol=2)
set.seed(998877)
alpha = .5
for (iboot in (1:nboot))
{
  un = ceiling((n-1)*runif(n-1))
  Rboot = R[un,]
  mean_vect = apply(Rboot,2,mean)
  grandMean = mean(mean_vect)
  mean_vect = alpha*mean_vect + (1-alpha)*grandMean   ## means shrunk toward grand mean
  cov_mat = cov(Rboot)
  sd_vect = sqrt(diag(cov_mat))
  Amat = cbind(rep(1,N),mean_vect,diag(1,N)) 
  muP = seq(min(mean_vect)+.001,max(mean_vect)-.001,length=300)                              
  sdP = muP 
  weights = matrix(0,nrow=300,ncol=N) 
  for (i in 1:length(muP))  
  {
    bvec = c(1,muP[i],rep(0,N))  
    result = 
      solve.QP(Dmat=2*cov_mat,dvec=rep(0,N),Amat=Amat,bvec=bvec,meq=2)
    sdP[i] = sqrt(result$value)
    weights[i,] = result$solution
  } 
  sharpe =( muP-mufree)/sdP 
  ind = (sharpe == max(sharpe)) 
  out[iboot,1] = sharpe[ind]
  wT = weights[ind,]
  sharpe_TRUE = (wT %*% mean_vect_TRUE - mufree) /
    sqrt(wT %*% cov_mat_TRUE %*% wT)
  out[iboot,2] = sharpe_TRUE
}
out_NoShortShrink = out
pdf("countries_boot_shrinkNoShrink.pdf",width=7,height=4)       ########  Figure 11.6  ########
par(mfrow=c(1,2))
boxplot(out_NoShort~gp,main="(a) No Shrinkage",ylim=c(0,.7))
abline(h=.3503,lwd=3,lty=2)
boxplot(out_NoShortShrink~gp,main="(b) Shrinkage",ylim=c(0,.7))
abline(h=.3503,lwd=3,lty=2)
graphics.off()

options(digits=3)
colMeans(out_NoShortShrink)
colMeans(out_NoShort)
sharpe_TRUE


################################################################
########## Code for figure 16.7  ###############################
################################################################

x = seq(0, 2, len = 300)
l1 = 0.25
l2 = 0.5
l3 = 1
l4 = 2
l5 = 5
c1 = 1 - exp(-l1)
c2 = 1 - exp(-l2)
c3 = 1 - exp(-l3)
c4 = 1 - exp(-l4)
c5 = 1 - exp(-l5)

pdf("utility_functions.pdf", width = 6, height = 5)
par(mfrow = c(1,1))
plot(x, (1 - exp(-l1*x))/c1, lwd = 2, type ="l", 
     ylab = expression(paste("U(x;",lambda,")") ) )
lines(x, (1 - exp(-l2*x))/c2, lwd = 2, col = "red", lty = 2)
lines(x, (1 - exp(-l3*x))/c3, lwd = 2, col = "blue", lty = 3)
lines(x, (1 - exp(-l4*x))/c4, lwd = 2, col = "purple", lty = 4)
lines(x, (1 - exp(-l5*x))/c5, lwd = 2, col = "cyan", lty = 4)
legend("topleft", c( expression(paste(lambda, "= 0.25")), 
                     expression(paste(lambda, "= 0.5")),
                     expression(paste(lambda, "= 1")),
                     expression(paste(lambda, "= 2")),
                     expression(paste(lambda, "= 5"))
), 
col=c("black", "red", "blue", "purple", "cyan"), lty = 1:5, lwd= 2   )
graphics.off()



################################################################
########## Code for Example 16.11  ################################
################################################################

library(quadprog)
dat = read.csv("Stock_Bond.csv")
y = dat[,c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21)]
n = dim(y)[1]
m = dim(y)[2] - 1
r = y[-1,]/y[-n,] - 1

mean_vect = as.matrix(colMeans(r))
cov_mat = cov(r)

nlambda = 250
loglambda_vect = seq(2, 8, length = nlambda)
w_matrix = matrix(nrow = nlambda, ncol = 10)
mu_vect = matrix(nrow = nlambda, ncol = 1)
sd_vect = mu_vect
ExUtil_vect = mu_vect
conv_vect = mu_vect
for (i in 1:nlambda)
{
  lambda = exp(loglambda_vect[i])
  opt = solve.QP(Dmat = as.matrix(lambda^2 * cov_mat), dvec = lambda * mean_vect, 
                 Amat = as.matrix(rep(1,10)), bvec = 1, meq = 1)
  w = opt$solution
  mu_vect[i] = w %*% mean_vect
  sd_vect[i] = sqrt(w %*% cov_mat %*% w)
  w_matrix[i,] = w
  ExUtil_vect[i] = opt$value
}

pdf("utility_frontier.pdf", width = 6.5, height = 2.75)  ###  Figure 16.8
par(mfrow = c(1, 3))
plot(loglambda_vect, mu_vect, type = "l", lwd = 2, 
     xlab = expression(paste("log(",lambda,")")), ylab = "E(return)"   )
plot(loglambda_vect, sd_vect, type = "l", lwd = 2, 
     xlab = expression(paste("log(",lambda,")")), ylab = "SD(return)"   )
plot(sd_vect, mu_vect, type = "l", lwd = 2, xlab = "SD(return)",
     ylab = "E(return)", main = "Efficient Frontier")
graphics.off()
