###########  Additional R functions for  #################################################
###########  Statistics and Data Analysis for Financial Engineering, 2nd Edition #########
###########  by Ruppert and Matteson  ####################################################

require(fGarch) 

#######################################
#### multivariate Ljung--Box tests ####
#######################################

"mq" = function(x, lag = 1, df.adj = 0){
  # Compute multivariate Ljung-Box test statistics
  #
  x = as.matrix(x)
  nr = dim(x)[1]
  nc = dim(x)[2]
  g0 = var(x)
  ginv = solve(g0)
  qm = 0.0
  df = 0
  out = matrix(0,lag,4)
  for (i in 1:lag){
    x1 = x[(i+1):nr,]
    x2 = x[1:(nr-i),]
    g = cov(x1,x2)
    g = g*(nr-i-1)/(nr-1)
    h = t(g)%*%ginv%*%g%*%ginv
    qm = qm+nr*nr*sum(diag(h))/(nr-i)
    df = df+nc*nc
    pv = 1-pchisq(qm,df-df.adj)
    #  print(c(i,qm,pv))
    out[i,] = c(i,round(qm,2),round(df-df.adj),round(pv,3))
  }
  output = as.data.frame(matrix(out[lag,],1,4))
  names(output) = c("K", "Q(K)", "d.f.", "p-value")                         
  print(output)
}

mLjungBox = mq

"mq0" = function(A, m){
  N = dim(A)[1]
  k = dim(A)[2]
  temp = numeric(m)
  pval = numeric(m)
  out = as.data.frame(matrix(0,m+1,4))
  names(out) = c("K", "Q(K)", "d.f.", "p-value")
  out[,1] = 0:m
  
  Q.temp = N*sum(cor(A)[lower.tri(cor(A), diag = FALSE)]^2)
  out[1,2] = Q.temp 
  df.temp = k*(k-1)/2
  out[1,3] = df.temp
  out[1,4] = 1-pchisq(Q.temp, df.temp)
  
  for(j in 1:m){
    ccf = cor(A[-(1:j),],A[-((N-j+1):N),])
    Q.temp = Q.temp + N*(N+2)*sum(ccf^2)/(N-j)
    out[(j+1),2] = Q.temp 
    df.temp = df.temp+ k^2
    out[(j+1),3] = df.temp
    out[(j+1),4] = 1-pchisq(Q.temp, df.temp)
  }  
  round(out,3)
}

###################################
#### EWMA estimation functions ####
###################################

"nllik.ewma" = function(lambda, innov){
  clambda = 1-lambda
  Sigma.hat = var(innov)
  invSigma.hat = chol2inv(chol(Sigma.hat)) # invSigma.hat = solve(Sigma.hat)
  detSigma.hat = det(Sigma.hat)
  llik = -0.5*log(detSigma.hat) - 0.5*crossprod(innov[1,],invSigma.hat)%*%innov[1,]
  llik = llik -0.5*log(detSigma.hat) - 0.5*crossprod(innov[2,],invSigma.hat)%*%innov[2,]
  n = dim(innov)[1]
  for(i in 3:n){
    atm1 = innov[(i-1),]
    at = innov[i,]
    denom = 1 - lambda^(i-1)
    #    Sigma.hat = clambda * tcrossprod(atm1) + lambda * Sigma.hat #approx
    Sigma.hat = (clambda/denom) * tcrossprod(atm1) + (lambda*(1-lambda^(i-2))/denom) * Sigma.hat #exact 
    invSigma.hat = chol2inv(chol(Sigma.hat)) # invSigma.hat = solve(Sigma.hat)
    detSigma.hat = det(Sigma.hat)
    llik = llik - 0.5*(log(detSigma.hat) + crossprod(at,invSigma.hat)%*%at)
  }
  nllik = -llik; nllik
}

"est.ewma" = function(lambda.0, innov){
  out = optim(lambda.0, nllik.ewma, lower = 0.001, upper = 0.999, innov = innov, method = "L-BFGS-B", hessian = TRUE)
  lambda.hat = out$par
  lambda.hat.se = 1/sqrt(out$hessian)
  list(lambda.hat = lambda.hat, lambda.hat.se = lambda.hat.se)
}

"sigma.ewma" = function(lambda, innov){
  clambda = 1-lambda
  Sigma.hat = var(innov)
  n = dim(innov)[1]
  d = dim(innov)[2]
  Sigma.t = array(0, c(d,d,n))
  Sigma.t[,,1:2] = Sigma.hat
  for(i in 3:n){
    atm1 = innov[(i-1),]
    at = innov[i,]
    denom = 1 - lambda^(i-1)
    #    Sigma.hat = clambda * tcrossprod(atm1) + lambda * Sigma.hat #approx
    Sigma.t[,,i] = (clambda/denom) * tcrossprod(atm1) + (lambda*(1-lambda^(i-2))/denom) * Sigma.t[,,(i-1)] #exact 
  }
  Sigma.t
}

###################################
#### DVEC GARCH(1,1) estimation ###
#### only for d = 2 components ####
#### requires library(mnormt) #####
###################################

"nllik.dvec.garch" = function(param, innov){
  omega = param[1:3]
  alpha = param[4:6]
  beta  = param[7:9]
  Y = as.matrix(innov)
  d = ncol(Y)
  n = nrow(Y)
  mu = numeric(d)
  Sigma.dvec = array(numeric(n*d*d),c(d,d,n))
  Sigma.dvec[,,1] = var(Y)
  V1 = Sigma.dvec[1,1,1]
  V2 = Sigma.dvec[2,2,1]
  V12 = Sigma.dvec[1,2,1]	
  llik = 0
  for(t in 2:n){
    V1  = omega[1] + alpha[1]*Y[(t-1),1]^2 + beta[1]*V1
    V2  = omega[2] + alpha[2]*Y[(t-1),2]^2 + beta[2]*V2
    V12 = omega[3] + alpha[3]*Y[(t-1),1]*Y[(t-1),2] + beta[3]*V12
    Sigma.dvec[1,1,t] = V1	
    Sigma.dvec[2,2,t] = V2
    Sigma.dvec[1,2,t] = V12
    Sigma.dvec[2,1,t] = V12
    llik = llik + dmnorm(Y[t,], mean = mu, varcov = Sigma.dvec[,,t], log = TRUE) 
  }
  -llik
}


"est.dvec.garch" = function(param, innov){
  out = optim(param, nllik.dvec.garch, lower = 0.0001, upper = c(rep(Inf,3),rep(0.9999,6)), innov = innov, method = "L-BFGS-B", hessian = TRUE,  control = list(trace = 6))
  theta.hat = out$par
  theta.hat.se = sqrt(diag(solve(out$hessian)))
  list(theta.hat = theta.hat, theta.hat.se = theta.hat.se)
}

"sigma.dvec.garch" = function(param, innov){
  omega = param[1:3]
  alpha = param[4:6]
  beta  = param[7:9]
  Y = as.matrix(innov)
  d = ncol(Y)
  n = nrow(Y)
  mu = numeric(d)
  Sigma.dvec = array(numeric(n*d*d),c(d,d,n))
  Sigma.dvec[,,1] = var(Y)
  V1 = Sigma.dvec[1,1,1]
  V2 = Sigma.dvec[2,2,1]
  V12 = Sigma.dvec[1,2,1]	
  for(t in 2:n){
    V1  = omega[1] + alpha[1]*Y[(t-1),1]^2 + beta[1]*V1
    V2  = omega[2] + alpha[2]*Y[(t-1),2]^2 + beta[2]*V2
    V12 = omega[3] + alpha[3]*Y[(t-1),1]*Y[(t-1),2] + beta[3]*V12
    Sigma.dvec[1,1,t] = V1	
    Sigma.dvec[2,2,t] = V2
    Sigma.dvec[1,2,t] = V12
    Sigma.dvec[2,1,t] = V12
  }
  list(Sigma.t = Sigma.dvec)
}

###################################
#### DCCe estimation functions ####
###################################

"llik.DCCe" = function(theta, innov){
  # llik for the correlation component
  E = innov
  n = dim(E)[1]
  d = dim(E)[2]
  
  lambda = theta[1]
  ####	alpha = theta[1]; beta = theta[2];	omega = (1-alpha-beta)
  
  S = cor(E)
  Q = S
  Q.temp = diag(d)
  
  llik = 0
  for(t in 2:n){
    Q = (1-lambda) * E[t-1,] %*% t(E[t-1,]) + lambda * Q
    ####		Q = S * omega + alpha * E[t-1,] %*% t(E[t-1,]) + beta * Q
    diag(Q.temp) = 1/sqrt(diag(Q))
    R = (Q.temp) %*% Q %*% (Q.temp)
    llik = llik + log(det(R)) + t(E[t,])%*%solve(R)%*%E[t,] - t(E[t,])%*%E[t,]
  }
  llik/2	
}


"est.DCCe" = function(theta, innov){
  out = optim(theta, llik.DCCe, lower = 0.001, upper = 0.999, innov = innov, method = "L-BFGS-B", hessian = FALSE)
  out$par
}

"sigma.DCCe" = function(theta, innov){
  Y = innov
  n = dim(Y)[1]
  d = dim(Y)[2]
  
  D = matrix(0,n,d)
  E = matrix(0,n,d)
  
  for(i in 1:d){
    fit.temp = garchFit(data = Y[,i], include.mean = FALSE, trace = FALSE)
    D[,i] = sqrt(fit.temp@h.t)
    E[,i] = Y[,i]/sqrt(fit.temp@h.t)
  }
  
  lambda = theta[1]
  ####	alpha = theta[1]; beta = theta[2];	omega = (1-alpha-beta)
  
  S = cor(E)
  Q = S
  Q.temp = diag(d)
  
  Sigma.t = array(0, c(d,d,n))
  R.t = array(0, c(d,d,n))
  Sigma.t[,,1] = var(Y)
  R.t[,,1] = S
  
  for(t in 2:n){
    Q = (1-lambda) * E[t-1,] %*% t(E[t-1,]) + lambda * Q
    ####		Q = S * omega + alpha * E[t-1,] %*% t(E[t-1,]) + beta * Q
    diag(Q.temp) = 1/sqrt(diag(Q))
    R.t[,,t] = (Q.temp) %*% Q %*% (Q.temp)
    Sigma.t[,,t] = diag(D[t,]) %*% R.t[,,t] %*% diag(D[t,])
  }
  
  list(Sigma.t = Sigma.t, R.t = R.t)
  
}

"fit.DCCe" = function(theta.0 = 0.9, innov){
  Y = innov
  n = dim(Y)[1]
  d = dim(Y)[2]
  
  D = matrix(0,n,d)
  E = matrix(0,n,d)
  
  garch.omega.alpha.beta = matrix(0,d,3)
  
  for(i in 1:d){
    fit.temp = garchFit(data = Y[,i], include.mean = FALSE, trace = FALSE)
    garch.omega.alpha.beta[i,] = fit.temp@fit$matcoef[,1]
    D[,i] = sqrt(fit.temp@h.t)
    E[,i] = Y[,i]/sqrt(fit.temp@h.t)
  }
  
  DCCe.params = est.DCCe(theta = theta.0, innov = E)
  DCCe.params
  DCCe.Sigma = sigma.DCCe(theta = DCCe.params, innov=Y)
  
  list(Sigma.t = DCCe.Sigma$Sigma.t, R.t = DCCe.Sigma$R.t, coef = garch.omega.alpha.beta, lambda = DCCe.params)
  
}

###################################
#### DOC estimation functions #####
###################################

"DOC.obj" = function(theta, Z, c, L){
  # Compute the unconstrained objective function for the DOC in volatility algotithm.
  # For use in esimating mixing parameter, theta (px1).
  # Time series vector Z (nxd) should be standardized.
  # Truncation level for computing Huber's function c (1x1) (non-negative).
  # Max_lag L (1x1) for lags 0,1,...,L. 
  # Diagonal weighting matrix.
  # See Matteson and Tsay (2010)
  p = length(theta)
  p2 = 2*p
  q = L*p2 + p
  N = L+1
  sizeZ = dim(Z)
  n = sizeZ[1]
  d = sizeZ[2]
  # Compute W from theta
  W = theta2W(theta)
  # Separate signals with estimate of W: S (nxd)
  S = Z %*% t(W)   
  # Transforming signals with Hubers function, saving column means
  SH = myHuber(S, c) 
  SH.bar = colMeans(SH)
  #Computing objective function
  f.bar = numeric(q)
  iter = 0:(2*L) * p + 1
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      SH.barXSH.bar = SH.bar[i] * SH.bar[j]
      f.bar[iter[1]]= crossprod(SH[,i],SH[,j])/n - SH.barXSH.bar
      for(ell in 1:L){
        f.bar[iter[2*ell]]   = crossprod(SH[-(1:ell),i],SH[-((n-ell+1):n),j])/n - SH.barXSH.bar
        f.bar[iter[2*ell+1]] = crossprod(SH[-(1:ell),j],SH[-((n-ell+1):n),i])/n - SH.barXSH.bar
      }
      iter = iter + 1
    }
  }  
  phi = 1 - (0:L)/N
  phi.total = sum(phi) 
  phi = phi / phi.total
  ##  PHI = c( rep(phi[1],p), rep(phi[-1],each = p2) )
  PHI = c( rep(phi[1]/p,p), rep(phi[-1]/p2,each = p2) )
  PHI = diag(PHI)
  crossprod(f.bar, PHI) %*% f.bar
}

"doc.garch" = function(E, L = 4., c = 2.25, theta.ini = NULL, 
                       n.ahead = 10, common.tail.index = FALSE){  
  E = as.matrix(E)
  d = ncol(E)
  n = nrow(E)
  p = d*(d-1)/2
  L = L 
  c = c 
  
  U = t(matrix.sqrt.inv(var(E)))
  U.inv = t(matrix.sqrt(var(E)))
  Z = E %*% U
  
  if(is.null(theta.ini)){
    theta.ini = as.matrix(rep(0,p)) 
  }
  else{
    theta.ini = as.matrix(theta.ini)
  }
  
  # Perform DOC in Volatility Estimation 
  out = optim(par= theta.ini, fn = DOC.obj, gr = NULL, 
              Z = Z, c = c, L = L, method = "L-BFGS-B",
              lower = -pi, upper = pi,
              control = list(), hessian = FALSE)
  W = canonicalW(theta2W(out$par))
  W.inv = t(W)
  theta.hat = W2theta(W)
  S = Z %*% t(W) 
  M.inv = U %*% t(W)
  M = W %*% U.inv
  
  # Esitmate sigma squared for DOCs and PCs:
  V.doc = array(numeric(n*d*d),c(d,d,n))
  V.pca = array(numeric(n*d*d),c(d,d,n))
  H.doc = array(numeric(n*d*d),c(d,d,n))
  H.pca = array(numeric(n*d*d),c(d,d,n))
  V.doc.forecast = array(numeric(n.ahead*d*d),c(d,d,n.ahead))
  V.pca.forecast = array(numeric(n.ahead*d*d),c(d,d,n.ahead))
  H.doc.forecast = array(numeric(n.ahead*d*d),c(d,d,n.ahead))
  H.pca.forecast = array(numeric(n.ahead*d*d),c(d,d,n.ahead))
  
  DOC.resid = matrix(0,n,d)
  PC.resid = matrix(0,n,d)
  
  DOC.est = matrix(0,d,4, dimnames = list(NULL,c("omega", "alpha1", "beta1", "shape")))
  PC.est = matrix(0,d,3, dimnames = list(NULL, 
                                         c("omega", "alpha1", "beta1")))
  
  if(common.tail.index == TRUE){
    tail.index.grid = seq(4.1,10,0.1)
    n.grid = length(tail.index.grid)
    tail.index.llh = numeric(n.grid)
    for(g in 1:n.grid){
      for(i in 1:d){
        fitDOC = garchFit(~garch(1,1), data=S[,i], 
                         include.mean = FALSE,
                         cond.dist =  "std", shape = tail.index.grid[g], include.shape = F,
                         trace=FALSE) #tail.index.grid[g], trace = F)
        tail.index.llh[g] = tail.index.llh[g] - fitDOC@fit$llh
     }
    }
    
    est.common.tail.index = tail.index.grid[which.max(tail.index.llh)]
    
    for(i in 1:d){
      fitDOC = garchFit(~garch(1,1), data=S[,i], 
                        include.mean = FALSE,
                        cond.dist =  "std", shape = est.common.tail.index,
                        include.shape = F, trace = F)
      
      DOC.resid[,i] = fitDOC@residuals / fitDOC@sigma.t
      V.doc[i,i,] = fitDOC@h.t
      V.doc.forecast[i,i,] = predict(fitDOC, n.ahead = n.ahead)$standardDeviation^2 
      DOC.est[i,] = c(fitDOC@fit$par, est.common.tail.index)
      
      fitPC = garchFit(~garch(1,1), data=Z[,i], 
                       include.mean = FALSE,
                       cond.dist =  "norm", trace = F)
      PC.resid[,i] = fitPC@residuals / fitPC@sigma.t               
      V.pca[i,i,] = fitPC@h.t
      V.pca.forecast[i,i,] = predict(fitPC, n.ahead = n.ahead)$standardDeviation^2 
      PC.est[i,] = fitPC@fit$par
    }
  }
 
  
  if(common.tail.index == FALSE){
   for(i in 1:d){
    fitDOC = garchFit(~garch(1,1), data=S[,i], 
                      include.mean = FALSE,
                      cond.dist =  "std", 
                      include.shape = TRUE, trace = F)
    
    DOC.resid[,i] = fitDOC@residuals / fitDOC@sigma.t
    V.doc[i,i,] = fitDOC@h.t
    V.doc.forecast[i,i,] = predict(fitDOC, n.ahead = n.ahead)$standardDeviation^2 
    DOC.est[i,] = fitDOC@fit$par
    
    fitPC = garchFit(~garch(1,1), data=Z[,i], 
                     include.mean = FALSE,
                     cond.dist =  "norm", trace = F)
    PC.resid[,i] = fitPC@residuals / fitPC@sigma.t               
    V.pca[i,i,] = fitPC@h.t
    V.pca.forecast[i,i,] = predict(fitPC, n.ahead = n.ahead)$standardDeviation^2 
    PC.est[i,] = fitPC@fit$par
   }
  }
  
  tM = t(M)
  tU.inv = t(U.inv)
  for(t in 1:n){
    H.doc[,,t] = tM %*% V.doc[,,t] %*% M
    H.pca[,,t] = tU.inv %*% V.pca[,,t] %*% U.inv
  }
  for(t in 1:n.ahead){
    H.doc.forecast[,,t] = tM %*% V.doc.forecast[,,t] %*% M
    H.pca.forecast[,,t] = tU.inv %*% V.pca.forecast[,,t] %*% U.inv   
  }
  #matrices W, U, M returned for dx1 (column) components s_t, z_t
  list(Sigma.doc = H.doc, Sigma.pca = H.pca, Z.hat = Z, S.hat = S, 
       coef.doc = DOC.est, coef.pca = PC.est, theta.hat = theta.hat,
       W.hat = W, U.hat = t(U), M.hat = tM, DOC.resid = DOC.resid, PC.resid = PC.resid,
       Sigma.doc.forecast = H.doc.forecast, Sigma.pca.forecast = H.pca.forecast)
  
} 

DOC.test = function(A, m){
  N = dim(A)[1]
  k = dim(A)[2]
  temp = numeric(m)
  pval = numeric(m)
  out = as.data.frame(matrix(0,m+1,4))
  names(out) = c("m", "Q(m)", "d.f.", "p-value")
  out[,1] = 0:m
  
  Q.temp = N*sum(cor(A)[lower.tri(cor(A), diag = FALSE)]^2)
  out[1,2] = Q.temp 
  df.temp = k*(k-1)/2
  out[1,3] = df.temp
  out[1,4] = 1-pchisq(Q.temp, df.temp)
  
  for(j in 1:m){
    ccf = cor(A[-(1:j),],A[-((N-j+1):N),])
    Q.temp = Q.temp + N*(N+2)*sum(ccf[lower.tri(ccf, diag = FALSE)]^2)/(N-j)
    Q.temp = Q.temp + N*(N+2)*sum(ccf[upper.tri(ccf, diag = FALSE)]^2)/(N-j)
    out[(j+1),2] = Q.temp 
    df.temp = df.temp+ k*(k-1)
    out[(j+1),3] = df.temp
    out[(j+1),4] = 1-pchisq(Q.temp, df.temp)
  }	
  round(out,3)
}


###################################
#### Other functions  #############
###################################

"matrix.sqrt" = function(A)
{
  sva <- svd(A)
  if (min(sva$d)>=0)
    Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
  else
    stop("Matrix square root is not defined")
  return(Asqrt)
}

"matrix.sqrt.inv" = function(A)
{
  sva <- svd(A)
  if (min(sva$d)>=0)
    Asqrt <- t(sva$v %*% (t(sva$u) / sqrt(sva$d)))
  else
    stop("Matrix square root is not defined")
  return(Asqrt)
}

"orthogonalize" <- function(M)
{
  solve(matrix.sqrt(M%*%t(M)))%*%M
}

pd.vcov.check = function(Sigma){
  dex = NULL
  if(is.na(dim(Sigma)[3])) {
    sva <- svd(Sigma)
    temp = min(sva$d)	
    if (min(sva$d)<=0) {dex = c(dex, t)}
  }
  else{ 
    N = dim(Sigma)[3] 
    temp = numeric(N)
    for(t in 1:N){
      sva <- svd(Sigma[,,t])
      temp[t] = min(sva$d)
      if (min(sva$d)<=0) {dex = c(dex, t)}
    }	
  }  		
  list(index = dex, min.sv = temp)
}



"mvwindow.var" <- function(rt,win=30){
  # rt: return series
  # win: window size.
  T=length(rt)
  vtemp = var(rt)
  vol=rep(vtemp,T)
  if (win < (T+1)){
    for (i in win:T){
      ist=i-win+1
      x=rt[ist:i]
      v=var(x)
      vol[i]=v
    }
  }
  list(volatility=vol)
}


"mvwindow.cor" <- function(A, B, win=30){
  # A,B: return series
  # win: window size.
  T=length(A)
  ctemp = cor(A,B)
  cor=rep(ctemp,T)
  if (win < (T+1)){
    for (i in win:T){
      ist=i-win+1
      x=A[ist:i]
      y=B[ist:i]
      v=cor(x,y)
      cor[i]=v
    }
  }
  list(correlation=cor)
}



"givens.rotation" <- function(theta=0, d=2, which=c(1,2))
{
  # For a given angle theta, returns a d x d Givens rotation matrix
  # Ex: for i < j , d = 2:  (c -s)
  #                         (s  c)
  c = cos(theta)
  s = sin(theta)
  M = diag(d)
  a = which[1]
  b = which[2]
  M[a,a] =  c
  M[b,b] =  c
  M[a,b] = -s
  M[b,a] =  s
  M
}

"theta2W" <- function(theta)
{
  # For a vector of angles theta, returns W, a d x d Givens rotation matrix:
  # W = Q.1,d %*% ... %*% Q.d-1,d %*% Q.1,d-1 %*% ... %*% Q.1,3 %*% Q.2,3 %*% Q.1,2 
  ##  if(theta < 0  || pi < theta){stop("theta must be in the interval [0,pi]")}
  d = (sqrt(8*length(theta)+1)+1)/2
  if(d - floor(d) != 0){stop("theta must have length: d(d-1)/2")}
  W = diag(d)
  index = 1
  for(j in 2:d){
    for(i in (j-1):1){
      Q.ij = givens.rotation(theta[index], d, c(i,j))
      W = Q.ij %*% W 
      index = index + 1
    }
  }
  W
}

"W2theta" <- function(W)
{
  # Decompose a d by d orthogonal matrix W into the product of
  # d(d-1)/2 Givens rotation matrices. Returns theta, the d(d-1)/2 by 1
  # vector of angles, theta.
  # W = Q.1,d %*% ... %*% Q.d-1,d %*% Q.1,d-1 %*% ... %*% Q.1,3 %*% Q.2,3 %*% Q.1,2 
  if(dim(W)[1] != dim(W)[2]){stop("W must be a square matrix")}
  d = dim(W)[1]
  #  if(sum(abs(t(W)%*%W  - diag(d))) > 1e-10){stop("W must be an orthogonal matrix")}
  theta = numeric(d*(d-1)/2)
  index = 1
  for(j in d:2){
    for(i in 1:(j-1)){      
      x = W[j,j]
      y = W[i,j]        
      theta.temp = atan2(y,x)    
      Q.ij = givens.rotation(theta.temp, d, c(i,j))
      W.temp = Q.ij %*% W  
      W = W.temp
      theta[index] = theta.temp
      index = index + 1
    }
  }
  -1*rev(theta)
}


"canonicalW" <- function(W)
{
  if(dim(W)[1] != dim(W)[2]){stop("W must be a square matrix")}
  d = dim(W)[1]
  #  if(sum(abs(t(W)%*%W  - diag(d))) > 1e-10){stop("W must be an orthogonal matrix")}
  W.temp = W
  W.new = matrix(0,d,d)
  for(i in 1:d){
    index = which.max(abs(W))
    row.index = index %% d #row
    row.index = ifelse(row.index == 0, d, row.index)
    col.index = ceiling(index / d) # col
    w.i = W.temp[row.index,]
    W.new[col.index,] = w.i * ifelse(w.i[col.index] < 0, -1, 1)
    W[row.index,] = 0
    W[,col.index] = 0
  }
  if(det(W.new) < 0) {W.new[d,] = -W.new[d,]}
  W.new
}


"myHuber" <- function(S, c){
  # Compute Huber's function with truncation level 'c'.
  if(c < 0){stop("truncation level 'c' must be non-negative")}
  S = abs(S)
  SH = S^2*ifelse(S <= c, 1,0) + (2*c*S - c^2)*ifelse(S  > c, 1,0)
  SH
}


###################################
#### DCCt estimation functions ####
###################################

"llik.DCCt" = function(theta, innov, m = 5){
  # llik for the correlation component
  E = innov
  n = dim(E)[1]
  d = dim(E)[2]
  
  theta1 = theta[1] 
  theta2 = theta[2]
  theta12 = (1-theta1-theta2)
  
  Gamma = cor(E)
  Gamma.t = Gamma
  
  llik = 0
  for(t in (m+1):n){
    Gamma.t = theta12 * Gamma + theta1 * Gamma.t + theta2 * cor(E[(t-1):(t-m),])
    llik = llik + log(det(Gamma.t)) + t(E[t,])%*%solve(Gamma.t)%*%E[t,] - t(E[t,])%*%E[t,]
  }
  llik/2	
}


"est.DCCt" = function(theta, innov, m = 5){
  out = optim(theta, llik.DCCt, lower = c(0.0001, 0.0001), upper = c(0.9999, 0.999), innov = innov, m = m, method = "L-BFGS-B", hessian = FALSE, control=list(trace=6))
  #    out = optim(theta, llik.DCCt, innov = innov, m = m, hessian = FALSE, control=list(trace=10))
  out$par
}

"sigma.DCCt" = function(theta, innov, m = 5){
  Y = innov
  n = dim(Y)[1]
  d = dim(Y)[2]
  
  D = matrix(0,n,d)
  E = matrix(0,n,d)
  
  for(i in 1:d){
    fit.temp = garchFit(data = Y[,i], include.mean = FALSE, trace = FALSE)
    D[,i] = sqrt(fit.temp@h.t)
    E[,i] = Y[,i]/sqrt(fit.temp@h.t)
  }
  
  theta1 = theta[1] 
  theta2 = theta[2]
  theta12 = (1-theta1-theta2)
  
  Gamma = cor(D)
  Gamma.t = Gamma
  
  Sigma.t = array(0, c(d,d,n))
  R.t = array(0, c(d,d,n))
  Sigma.t[,,1:m] = var(Y)
  R.t[,,1:m] = Gamma
  
  for(t in (m+1):n){
    Gamma.t = theta12 * Gamma + theta1 * Gamma.t + theta2 * cor(D[(t-1):(t-m),])
    R.t[,,t] = Gamma.t
    Sigma.t[,,t] = diag(D[t,]) %*% R.t[,,t] %*% diag(D[t,])
  }
  
  list(Sigma.t = Sigma.t, R.t = R.t)
  
}

"fit.DCCt" = function(theta.0 = 0.9, innov, m = 5){
  Y = innov
  n = dim(Y)[1]
  d = dim(Y)[2]
  
  D = matrix(0,n,d)
  E = matrix(0,n,d)
  
  garch.omega.alpha.beta = matrix(0,d,3)
  
  for(i in 1:d){
    fit.temp = garchFit(data = Y[,i], include.mean = FALSE, trace = FALSE)
    garch.omega.alpha.beta[i,] = fit.temp@fit$matcoef[,1]
    D[,i] = sqrt(fit.temp@h.t)
    E[,i] = Y[,i]/sqrt(fit.temp@h.t)
  }
  
  DCCt.params = est.DCCt(theta = theta.0, innov = E, m = m)
  DCCt.params
  DCCt.Sigma = sigma.DCCt(theta = DCCt.params, innov=Y, m = m)
  
  list(Sigma.t = DCCt.Sigma$Sigma.t, R.t = DCCt.Sigma$R.t, coef = garch.omega.alpha.beta, lambda = DCCt.params)
  
}

