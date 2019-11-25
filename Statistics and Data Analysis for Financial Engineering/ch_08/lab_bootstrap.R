library(bootstrap)
library(MASS)

set.seed("3857")

data(CRSPday, package = "Ecdat")

ge <- CRSPday[, 4]

getRetVsT <- function(returns, df_canidates = c(1, 2, 4, 6, 10, 20)) {
  plots <- lapply(df_canidates, function(df) {
    
    n <- length(returns)
    q_range <- (1:n) / (n+1)
    
    data <- data.table(ret = returns, theoretical = qt(q_range, df))
    data$theoretical <- sort(data$theoretical)
    
    model <- lm(qt(c(0.25,0.75), df = df) ~ quantile(data$ret,c(0.25,0.75)))
    
    ggplot(data, aes(x = sort(ret), y = theoretical)) +
      geom_abline(col = 'cornflowerblue', lwd = 1.3, slope = model$coefficients[2], intercept = model$coefficients[1]) +
      geom_point() +
      labs(title = paste("df = ", df), 
           x = "returns",
           y = "theoretical")
  })
  
  do.call(grid.arrange, c(plots, top = "QQ-Plot: returns vs t-distribution"))
}

getRetVsT(ge)

nboot <- 1000

t_mle <- function(x) { as.vector(
  fitdistr(x, "t", start = list( m = mean(x), s = sd(x), df = 3) )$estimate
)}

results <- bootstrap(ge, nboot, t_mle)

rowMeans(results$thetastar[ , ])
apply(results$thetastar[,], 1, sd)

fitdistr(ge, "t")
