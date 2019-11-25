library(copula)

u <- seq(0.000001, 1, length = 500)
frank <- iPsi(copula = archmCopula(family = "frank", param = 1), u)
plot(u, frank, type = "l", lwd = 3, ylab = expression(phi(u)))
abline(h = 0); abline(v = 0)

set.seed(5640)

theta <- c(-100, -50, -10, -1, 0, 5, 20, 50, 500)
par(mfrow = c(3,3), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)

for(i in 1:9) {
  U <- rCopula(n = 200,
               copula = archmCopula(family = "frank", param = theta[i]))
  plot(U, xlab = expression(u[1]), ylab = expression(u[2]),
       main = eval(substitute(expression(paste(theta, " = ", j)),
                              list(j = as.character(theta[i])))))
}

set.seed(5640)

theta <- c(-0.98, -0.7, -0.3, -0.1, 0.1, 1, 1.5, 15, 100)
par(mfrow = c(3,3), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)

for(i in 1:9) {
  U <- rCopula(n = 200,
               copula = archmCopula(family = "clayton", param = theta[i]))
  plot(U, xlab = expression(u[1]), ylab = expression(u[2]),
       main = eval(substitute(expression(paste(theta, " = ", j)),
                              list(j = as.character(theta[i])))))
}

set.seed(5640)

theta <- c(1.1, 1.5, 2, 4, 8, 50)
par(mfrow = c(2,3), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)

for(i in 1:6) {
  U <- rCopula(n = 200,
               copula = archmCopula(family = "gumbel", param = theta[i]))
  plot(U, xlab = expression(u[1]), ylab = expression(u[2]),
       main = eval(substitute(expression(paste(theta, " = ", j)),
                              list(j = as.character(theta[i])))))
}

set.seed(5640)

theta <- c(1.1, 1.5, 2, 4, 8, 50)
par(mfrow = c(2,3), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)

for(i in 1:6) {
  U <- rCopula(n = 200,
               copula = archmCopula(family = "gumbel", param = theta[i]))
  plot(U, xlab = expression(u[1]), ylab = expression(u[2]),
       main = eval(substitute(expression(paste(theta, " = ", j)),
                              list(j = as.character(theta[i])))))
}

rho <- seq(-1, 1, by = 0.01)
df <- c(1, 4, 25, 240)

x1 <- -sqrt((df[1] + 1) * (1-rho)/(1+rho))
lambda1 <- 2 * pt(x1, df[1] + 1)

x4 <- -sqrt((df[2]+1)*(1-rho)/(1+rho))
lambda4 <- 2 * pt(x4, df[2] + 1)

x25 <- -sqrt((df[3]+1)*(1-rho)/(1+rho))
lambda25 <- 2 * pt(x25, df[3] + 1)

x250 <- -sqrt((df[4]+1)*(1-rho)/(1+rho))
lambda250 <- 2 * pt(x250, df[4] + 1)

par(mfrow = c(1,1), lwd = 2, cex.axis = 1.2, cex.lab = 1.2)
plot(rho, lambda1, type = "l", lty = 1, xlab = expression(rho),
     ylab = expression(lambda[l] == lambda[u]))
lines(rho, lambda4, lty=2)
lines(rho, lambda25, lty = 3)
lines(rho, lambda250, lty = 4)
legend("topleft", c(expression(nu==1), expression(nu==4),
                    expression(nu==25), expression(nu == 250)), lty = 1:4)

