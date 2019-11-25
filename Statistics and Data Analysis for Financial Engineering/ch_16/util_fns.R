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
# Computes the Efficient Frontier.
#
#-----

library(quadprog)
library(linprog)

efficient_frontier = function(returns, muP, mufree=0.0, w_lower_limit=-Inf, w_upper_limit=+Inf){
  #
  # w_lower_limit and w_upper_limit (if finite) are bounds on the porfolio weights such that
  # 
  # w_lower_limit <= w <= w_upper_limit
  #
  # Note: if solve.QP gives errors about convergence one can try a more restricted range for the values of muP:
  #
  # muP = seq(min(mean_vect), max(mean_vect), length.out=500)
  #
  #--
  n_stocks = dim(returns)[2]
  
  # Extract individual equity (mean return, std return) values: 
  mean_vect = apply(returns, 2, mean)
  cov_mat = cov(returns)
  sd_vect = sqrt(diag(cov_mat))
  
  # Portfolio weight constraints:
  #
  #     first condition enforce \sum w_j = 1
  #     second condition enforce \sum mu_j w_j = mu_target
  #             
  Amat = cbind(rep(1,n_stocks),mean_vect)
  bvec = c(1,NaN)
  if( is.finite(w_lower_limit) ){
    Amat = cbind(Amat, diag(1, nrow=n_stocks))
    bvec = c(bvec, w_lower_limit*rep(1,n_stocks))
  }
  if( is.finite(w_upper_limit) ){
    Amat = cbind(Amat, -diag(1, nrow=n_stocks))
    bvec = c(bvec, -w_upper_limit*rep(1,n_stocks))
  }

  # storage for results: 
  sdP = muP
  weights = matrix(0, nrow=length(muP), ncol=n_stocks)
  
  print(bvec)
  print(Amat)
  
  # find the optimal portfolios for each target expected return:
  #
  for( i in 1:length(muP) ){
    bvec[2] = +muP[i] # enforce portfolio mean constraint 
    result = solve.QP(Dmat=2*cov_mat, dvec=rep(0,n_stocks), Amat=Amat, bvec=bvec, meq=2)
    sdP[i] = sqrt(result$value)
    weights[i,] = result$solution 
  }

  # Find maximum Sharpe portfolio: 
  sharpe = ( muP - mufree ) / sdP 
  ind_ms = which.max(sharpe)

  # Find minimum variance portfolio: 
  ind_mv = which.min(sdP)
  
  list(mufree=mufree, mean_vect=mean_vect, cov_mat=cov_mat, sd_vect=sd_vect,
       muP=muP, sdP=sdP, weights=weights,
       sharpe=sharpe, max_sharpe=sharpe[ind_ms], max_sharpe_ind=ind_ms,
       min_variance=sdP[ind_mv]^2, min_variance_ind=ind_mv )
}


plot_efficient_frontier = function(result, title=''){
    #
    # Plots the output from a call to the function "efficient_frontier".
    #
    # weights[ind_ms, ] # print the weights of the tangency portfolio
    #

    # Unpack:
    #
    mufree = result$mufree
    muP = result$muP
    sdP = result$sdP
    weights = result$weights
    sharpe = result$sharpe
    ind_ms = result$max_sharpe_ind
    ind_mv = result$min_variance_ind
    mean_vect = result$mean_vect
    sd_vect = result$sd_vect
    
    max_sdP = max(c(sdP, sd_vect))
    x_lim_max = 1.05*max_sdP
    max_muP = max(c(muP, mean_vect))
    y_lim_max = 1.05*max_muP
    
    plot(sdP, muP, type='l', xlim=c(0, x_lim_max), ylim=c(0, y_lim_max), xlab='sdP', ylab='muP', main=title)
    ind3 = (muP > muP[ind_mv])
    lines(sdP[ind3], muP[ind3], type='l', lwd=3, col='red')
    points(0, mufree, cex=4, pch='*')
    grid()

    # Plot the tangency portfolio:
    #
    sdP_max = sdP[ind_ms]
    lines(c(0, sdP_max), mufree + c(0, sdP_max) * (muP[ind_ms] - mufree)/sdP[ind_ms], lwd=4, lty=1, col='blue')
    points(sdP[ind_ms], muP[ind_ms], cex=4, pch='*') # tangency portfolio

    points(sdP[ind_mv], muP[ind_mv], cex=2, pch='+') # min variance portfolio
}


maximize_expected_utility = function(returns, loglambda_vect, mufree=0.0, w_lower_limit=-Inf, w_upper_limit=+Inf){
    #
    # w_lower_limit and w_upper_limit (if finite) are bounds on the porfolio weights such that
    #
    # w_lower_limit <= w <= w_upper_limit
    #
    # Note: if solve.QP gives errors about convergence one can try a more restricted range for the values of muP:
    #
    # muP = seq(min(mean_vect), max(mean_vect), length. out=500)
    #
    #---
    n_stocks = dim(returns)[2]
    
    # Extract individual equity (mean return, std return) values:
    mean_vect = apply(returns, 2, mean)
    cov_mat = cov(returns)
    sd_vect = sqrt(diag(cov_mat))

    # Portfolio weight constraints:
    #
    # first condition enforces \sum w_j = 1
    #
    Amat = matrix(rep(1, n_stocks), nrow=n_stocks, ncol=1)
    bvec = c(1)
    if( is.finite(w_lower_limit) ){
        Amat = cbind(Amat, diag(1, nrow=n_stocks) )
        bvec = c(bvec, w_lower_limit*rep(1,n_stocks) )
    }
    if( is.finite(w_upper_limit) ){
        Amat = cbind(Amat, -diag(1, nrow=n_stocks) )
        bvec = c(bvec, -w_upper_limit*rep(1,n_stocks) )
    }
    
    # storage for results:
    weights = matrix(0, nrow=length(loglambda_vect), ncol=n_stocks)
    mu_vect = matrix(0, nrow=length(loglambda_vect), ncol=1)
    sd_vect = matrix(0, nrow=length(loglambda_vect), ncol=1)
    ExUtil_vect = matrix(0, nrow=length(loglambda_vect), ncol=1)
    
    # find the optimal portfolios for each value of lambda:
    #
    for( i in 1:length(loglambda_vect) ){
        lambda = exp(loglambda_vect[i])
        result = solve.QP(Dmat=lambda*cov_mat, dvec=mean_vect, Amat=Amat, bvec=bvec, meq=1)
        w = result$solution
        mu_vect[i] = w %*% mean_vect
        sd_vect[i] = sqrt(w %*% cov_mat %*% w)
        weights[i,] = w
        ExUtil_vect[i] = result$value
    }
    
    list(mean_vect=mean_vect, cov_mat=cov_mat, sd_vect=sd_vect,
         weights=weights, mu_vect=mu_vect, sd_vect=sd_vect, ExUtil_vect=ExUtil_vect)

}


possible_expected_returns = function(returns, B1, B2){
    #
    # Note: B1 is the UPPER bound on the weights in the portfolio and
    #      -B2 is the LOWER bound on the weights in the portfolio that is:
    #
    # -B2 <= w_i <= +B1
    #
    stopifnot (B1>=0)
    stopifnot (B2>=0)

    mean_vect = colMeans(returns)
    M = length(mean_vect)

    AmatLP1 = cbind(diag(1, nrow=M), matrix(0, nrow=M, ncol=M))
    AmatLP2 = cbind(matrix(0, nrow=M, ncol=M), diag(1, nrow=M))
    AmatLP3 = c(rep(1, M), rep(-1, M))
    AmatLP = rbind(AmatLP1, AmatLP2, AmatLP3)
    bvecLP = c(rep(B1, M), rep(B2, M), 1)
    cLP = c(mean_vect, -mean_vect)
    const.dir = c(rep('<=', 2*M), '=')

    resultLP_min = solveLP(cvec=cLP, bvec=bvecLP, Amat=AmatLP, lpSolve=TRUE, const.dir=const.dir, maximum=FALSE)
    resultLP_max = solveLP(cvec=cLP, bvec=bvecLP, Amat=AmatLP, lpSolve=TRUE, const.dir=const.dir, maximum=TRUE)

    list(minRet=resultLP_min, maxRet=resultLP_max)
}

