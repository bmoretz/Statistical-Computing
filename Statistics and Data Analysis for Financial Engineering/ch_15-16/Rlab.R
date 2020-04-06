#
# Written by:
# -- 
# John L. Weatherwax                2009-04-21
# 
# email: wax@alum.mit.edu
# 
# Please send comments and especially bug reports to the
# above email address.
# 
#-----

if( !require('linprog') ){
    install.packages('linprog', dependencies=TRUE, repos='http://cran.rstudio.com/')
}
library(linprog)

save_plots = FALSE

set.seed(0)

src.dir <- "D:/Projects/MSDS-RiskAnalytics/Module_08/"

setwd(src.dir)

source('util_fns.R')

# Test the function "efficient_frontier" using the book's example on EPage 478:
#
library(Ecdat)
data(CRSPday)
R = 100*CRSPday[,4:6] # in percents
mean_vect = apply(R, 2, mean)
cov_mat = cov(R)
sd_vect = sqrt(diag(cov_mat))

# The target portfolio returns:
muP = seq(0.05, 0.14, length.out=300)
mufree = 1.3/253 

result = efficient_frontier(R, muP, mufree)

mean_vect = result$mean_vect
sd_vect = result$sd_vect
weights = result$weights
ind_ms = result$max_sharpe_ind

# Plot our results:
#
plot_efficient_frontier(result, "Duplicate Fig. 16.3")

print(weights[ind_ms, ]) # print the weights of the tangency portfolio

text(sd_vect[1], mean_vect[1], 'GE', cex=1.5)
text(sd_vect[2], mean_vect[2], 'IBM', cex=1.5)
text(sd_vect[3], mean_vect[3], 'Mobil', cex=1.5)



# Test the function "maximize_expected_utility" using the book's example on EPage 487:
#
dat = read.csv("../../BookCode/Data/Stock_Bond.csv", header=TRUE)
y = dat[, c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21)]
n = dim(y)[1]
m = dim(y)[2]
r= y[-1,] / y[-n, ] - 1

nlambda = 250
loglambda_vect = seq(2, 8, length=nlambda)

meu_result = maximize_expected_utility(r, loglambda_vect)


# This code reproduces the output of Fig. 16.8:
#
par(mfrow=c(1, 3))
#
plot(loglambda_vect, meu_result$mu_vect, type='l', xlab='log(\u03BB)', ylab='E(return)') # \u03BB is the unicode character for the LaTex symbol $\lambda$
grid()
plot(loglambda_vect, meu_result$sd_vect, type='l', xlab='log(\u03BB)', ylab='SD(return)')
grid()
plot(meu_result$sd_vect, meu_result$mu_vect, type='l', xlab='SD(return)', ylab='E(return)', main='Efficient Frontier')
grid()
par(mfrow=c(1, 1))


# P 1: EPage 488-489:
#
dat = read.csv("../../BookCode/Data/Stock_Bond.csv", header=TRUE)
prices = cbind(dat$GM_AC, dat$F_AC, dat$CAT_AC, dat$UTX_AC, dat$MRK_AC, dat$IBM_AC)
n = dim(prices)[1]
returns = 100 * ( prices[2:n,] / prices[1:(n-1),] - 1 )
#pairs(returns)

mean_vect = colMeans(returns)
cov_mat = cov(returns)
sd_vect = sqrt(diag(cov_mat) )

mufree = 3/365

muP = seq(min(mean_vect), max(mean_vect), length.out=500)

result = efficient_frontier(returns, muP, mufree, w_lower_limit=-0.1, w_upper_limit=0.5)

mean_vect = result$mean_vect
sd_vect = result$sd_vect
weights = result$weights
ind_ms = result$max_sharpe_ind

if( save_plots ){ postscript("../../WriteUp/Graphics/Chapter16/chap_16_rlab_prob_1_ef.eps", onefile=FALSE, horizontal=FALSE) }

plot_efficient_frontier(result, 'Problem 1')

print(round(weights[ind_ms,], 4)) # print the weights of the tangency portfolio

text(sd_vect[1], mean_vect[1], 'GM', cex=1.5)
text(sd_vect[2], mean_vect[2], 'F', cex=1.5)
text(sd_vect[3], mean_vect[3], 'CAT', cex=1.5)
text(sd_vect[4], mean_vect[4], 'UTX', cex=1.5)
text(sd_vect[5], mean_vect[5], 'MRK', cex=1.5)
text(sd_vect[6], mean_vect[6], 'TBM', cex=1.5)

if( save_plots ){ dev.off() }


# P 2: EPage 489
#
ind = result$max_sharpe_ind

E_R_P = 0.0007 * 100 # the desired return (as a percent)
E_R_T = result$muP[ind] # the return of the tangency portfolio
c( E_R_P, E_R_T, mufree )

omega = ( E_R_P - mufree ) / ( E_R_T - mufree )
print(sprintf("omega= %10.6f", omega) )

print("Tangency porfolio weights:")
ind = result$max_sharpe_ind
print(round(omega*result$weights[ind,], 4))

S = 100000
S * omega * result$weights[ind, ]


# P 4; EPage 489-490
#
dat = read.csv("../../BookCode/Data/FourStocks_Daily2013.csv", header=TRUE)
prices = dat[, -1]
n = dim(prices)[1]
returns = 100 * ( prices[-1,] / prices[-n,] - 1 )

mufree = 1.3/365

muP = seq(0.045, 0.06, length.out=500)

result = efficient_frontier(returns, muP, mufree, w_lower_limit=-0.5, w_upper_limit=0.5)

mean_vect = result$mean_vect
sd_vect = result$sd_vect
weights = result$weights
ind_ms = result$max_sharpe_ind

if( save_plots ){ postscript("../../WriteUp/Graphics/Chapter16/chap_16_rlab_prob_4_ef.eps", onefile=FALSE, horizontal=FALSE) }

plot_efficient_frontier(result, 'Problem 4')

print(sprintf('max sharpe (daily)= %10.6f', result$max_sharpe) )
print(sprintf('max sharpe (yearly)= %10.6f', result$max_sharpe * sqrt(252)))
print('Tangency portfolio weights:')
print(round(weights[ind_ms,], 4)) # print the weights of the tangency portfolio

text(sd_vect[1], mean_vect[1], 'AAPL', cex=1.5)
text(sd_vect[2], mean_vect[2], 'XOM', cex=1.5)
text(sd_vect[3], mean_vect[3], 'TGT', cex=1.5)
text(sd_vect[4], mean_vect[4], 'MCD', cex=1.5)

if( save_plots ){ dev.off() }


# P 5: EPage 490
#
dat = read.csv("../../BookCode/Data/FourStocks_Daily2013.csv", header=TRUE)
prices = dat[, -1]
n = dim(prices)[1]
returns = 100 * ( prices[-1,] / prices[-n,] - 1 )

loglambda_vect = seq(1.e-3, 8, length.out=200)

meu_result = maximize_expected_utility(returns, loglambda_vect)

if( save_plots ){
    postscript("../../WriteUp/Graphics/Chapter16/chap_16_rlab_prob_5.eps", onefile=FALSE, horizontal=FALSE)
    xlab='log(lambda)'
}else{
    xlab='log(\u03BB)' # \u03BB is the unicode character for the LaTex symbol $\lambda$
}
par(mfrow=c(1, 3))
plot(loglambda_vect, meu_result$mu_vect, type='l', xlab=xlab, ylab='E(return)') 
abline(h=0.0506, col='red')
grid()
plot(loglambda_vect, meu_result$sd_vect, type='l', xlab=xlab, ylab='SD(return)')
grid()
plot(meu_result$sd_vect, meu_result$mu_vect, type='l', xlab='SD(return)', ylab='E(return)', main='Efficient Frontier')
grid()
par(mfrow=c(1, 1))
if( save_plots ){ dev.off() }


mu_P = 0.0506
indx = which.min(abs(meu_result$mu_vect - mu_P))
print(sprintf('For mu_P= %f take lambda= %f', mu_P, exp(loglambda_vect[indx])))

#indx = which.min(abs(meu_result$sd_vect - sqrt(result$min_variance) ) )
#print(sprintf('closest I can get to min_sd= %f is sf where lambda= %f', sqrt(result$min_variance), sqrt(meu_result$sd_vect[indx]), exp(loglambda_vect[indx]) ))


dat = read.csv("../../BookCode/Data/Stock_Bond.csv", header=TRUE)
prices = cbind(dat$GM_AC, dat$F_AC, dat$CAT_AC, dat$UTX_AC, dat$MRK_AC, dat$IBM_AC)
n = dim(prices)[1]
returns = 100 * ( prices[2:n,] / prices[1:(n-1),] - 1 )


# As a "test" of the routine possible_expected_returns lets verify that a long only portfolio has its min/max returns
# given by the smallest/largest returns from the assets
#
possible_rets = possible_expected_returns(returns, B1=1.0, B2=0.0)
print('colMeans(returns)= ')
print(colMeans(returns))
print(sprintf('Long only portfolio return bounds: minRet= %f; maxRet= %f', possible_rets$minRet$opt, possible_rets$maxRet$opt))


# P 6: EPage 490-491
#
possible_rets = possible_expected_returns(returns, B1=0.3, B2=0.1)
print(sprintf('Problem 6: minRet= %f; maxRet= %f', possible_rets$minRet$opt, possible_rets$maxRet$opt))


# P 7: EPage 491
#
possible_rets = possible_expected_returns(returns, B1=0.15, B2=0.15)
print(possible_rets) # note both solveLP report errors
