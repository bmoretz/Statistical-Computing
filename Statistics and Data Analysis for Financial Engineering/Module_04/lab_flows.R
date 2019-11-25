library(data.table)
library(ggplot)
library(ggtheme)
library(Ecdat)
library(faraway)
library(fGarch)
library(sn)
library(copula)

theme_set(theme_light())

setwd("D:/Projects/MSDS-RiskAnalytics/datasets")

# Theme Overrides
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             legend.position = "top", legend.title = element_blank())

dat <- read.csv("FlowData.csv")
dat <- dat / 1000

n <- nrow(dat)

x1 <- dat$Flow1
fit1 <- st.mple(matrix(1, n, 1), y = x1, dp = c(mean(x1), sd(x1), 0, 10))
est1 <- fit1$dp
u1 <- pst(x1, dp = est1)

x2 <- dat$Flow2
fit2 <- st.mple(matrix(1, n, 1), y = x2, dp = c(mean(x2), sd(x2), 0, 10))
est2 <- fit2$dp

u2 <- pst(x2, dp = est2)

U.hat <- cbind(u1, u2)

z1 <- qnorm(u1)
z2 <- qnorm(u2)

Z.hat <- cbind(z1, z2)

library(ks)

fhatU <- kde(x = U.hat, H = Hscv(x=U.hat))
par(mfrow = c(2, 2), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)
hist(u1, main = "(a)", xlab = expression(hat(U)[1]), freq = F)
hist(u2, main = "(b)", xlab = expression(hat(U)[2]), freq = F)
plot(u1, u2, main = "(c)", xlab = expression(hat(U)[1]),
     ylab = expression(hat(U)[2]), mgp = c(2.5, 1, 0))
plot(fhatU, drawpoints = F, drawlabels = F,
     cont = seq(10, 80, 10), main = "(d)", xlab = expression(hat(U)[1]),
     ylab = expression(hat(U)[2]), mgp = c(2.5, 1, 0))


cor.test(u1, u2, method = "spearman")
cor.test(u1, u2, method = "kendall")
sin(-0.242*pi/2)
cor.test(u1, u2, method = "pearson")
cor.test(z1, z2, method = "pearson")

omega <- -0.371

#
Ct <- fitCopula(copula = tCopula(dim = 2), data = U.hat,
                method = "ml", start = c(omega, 10))
Ct@estimate

loglikCopula(param = Ct@estimate, u = U.hat, copula = tCopula(dim = 2))
-2*.Last.value + 2*length(Ct@estimate)
#
Cgauss <- fitCopula(copula = normalCopula(dim = 2), data = U.hat,
                    method = "ml", start = c(omega))
Cgauss@estimate

loglikCopula(param = Cgauss@estimate, u = U.hat,
             copula = normalCopula(dim = 2))
-2*.Last.value + 2*length(Cgauss@estimate)
#
Cfr <- fitCopula(copula = frankCopula(1, dim = 2), data = U.hat,
                 method = "ml")
Cfr@estimate
loglikCopula(param = Cfr@estimate, u = U.hat,
             copula = frankCopula(dim = 2))
-2*.Last.value + 2 * length(Cfr@estimate)
#
Ccl <- fitCopula(copula = claytonCopula(1, dim = 2), data = U.hat,
                 method = "ml")
Ccl@estimate
loglikCopula(param = Ccl@estimate, u = U.hat,
             copula = claytonCopula(dim = 2))
-2*.Last.value + 2*length(Ccl@estimate)

pdf("unif_flows_contours_copulas.pdf",width=7,height=6)
#
par(mfrow=c(2,3), mgp = c(2.5, 1, 0))
plot(u1, u2, main="Uniform-Transformed Data",
     xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))
Udex = (1:n)/(n+1)
Cn = C.n(u = cbind(rep(Udex, n), rep(Udex, each=n)) , X = U.hat, 
         offset=0)
EmpCop = expression(contour(Udex,Udex,matrix(Cn,n,n), col=2, add=T))
#
contour(normalCopula(param=0,dim=2), pCopula, main=expression(C[0]), 
        xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))
eval(EmpCop)
#
contour(tCopula(param=Ct@estimate[1], dim=2, 
                df=round(Ct@estimate[2])), 
        pCopula, main = expression(hat(C)[t]), 
        xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))
eval(EmpCop)
#
contour(normalCopula(param=Cgauss@estimate[1], dim = 2), 
        pCopula, main = expression(hat(C)[Gauss]),
        xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))
eval(EmpCop)
#
contour(frankCopula(param=Cfr@estimate[1], dim = 2), 
        pCopula, main = expression(hat(C)[Fr]),
        xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))
eval(EmpCop)
#
contour(claytonCopula(param=Ccl@estimate[1], dim = 2), 
        pCopula, main = expression(hat(C)[Cl]),
        xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))
eval(EmpCop)
#
graphics.off()
