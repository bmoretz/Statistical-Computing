library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(GGally)
library(ggthemes)
library(scales)
library(reshape2)
library(skimr)
library(gridExtra)
library(Ecdat)
library(faraway)
library(fGarch)

#####################################################################
######################### Ford Case Study ###########################
#####################################################################

theme_set(theme_sjplot())

# Theme Overrides
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             legend.position = "top", legend.title = element_blank())

path.data <- "D:/Projects/MSDS-RiskAnalytics/datasets"
setwd(path.data)

# 4.11
# Ford Case Study

data.ford <- as.data.table(read.csv("ford.csv", header = T))

colnames(data.ford) <- c("n", "date", "return")

head(data.ford)

# E1

mean(data.ford$return)
sd(data.ford$return)
median(data.ford$return)

ggplot(data.ford, aes(sample = return)) +
  geom_qq() +
  geom_qq_line()

shapiro.test(data.ford$return)

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

getRetVsT(data.ford$return)

F_inv_q <- median(data.ford$return)
d <- density(data.ford$return)
a <- approx(d$x, d$y, F_inv_q, method = "linear")
q <- 0.5

sqrt( ( q*(1-q) ) / (length(data.ford$return) * a$y^2) )

sd(data.ford$return) / sqrt(length(data.ford$return))

# E2

getRetDensityVsNorm <- function(returns) {
  data <- data.table(x = seq(min(returns), max(returns), length.out = length(returns)), y = returns)
  
  ggplot(data, aes(y)) +
    geom_density(aes(col = "KDE"), lwd = 1) +
    geom_line(aes(x, dnorm(x, mean = mean(returns), sd = sd(returns)), col = "normal(mean, sd)"), lwd = 1.1, linetype = "longdash") +
    geom_line(aes(x, dnorm(x, mean = median(returns), sd = mad(returns)), col = "normal(median, mad)"), lwd = 1.1, linetype = "dashed") +
    scale_y_continuous(labels = scales::comma) +
    labs(title = paste("KDE: Returns Vs Normal, n=", nrow(data)), y = "density", x = "") +
    theme(legend.position = "right")
}

getRetDensityVsNorm(diff(Garch$dy))

getRetNormQuantiles <- function(returns, quantiles = c(0.25,0.1,0.05,0.025,0.01,0.0025), desc = "") {
  plots <- lapply(quantiles, function(p) {
    
    p_value <- p
    n <- length(returns)
    q_range <- (1:n) / (n+1)
    
    data <- data.table(ret = returns, theoretical = qnorm(q_range))
    data$theoretical <- sort(data$theoretical)
    
    model <- lm(qnorm(c(p_value, 1-p_value)) ~ quantile(data$ret, c(p_value, 1-p_value)))
    
    ggplot(data, aes(x = sort(ret), y = theoretical)) +
      geom_abline(col = 'cornflowerblue', lwd = 1.3, slope = model$coefficients[2], intercept = model$coefficients[1]) +
      geom_point() +
      labs(title = paste("p = ", p_value), 
           x = "returns",
           y = "theoretical")
  })
  
  do.call(grid.arrange, c(plots, top = paste(desc, "vs Normal Quantiles")))
}

getRetNormQuantiles(diff(Garch$bp), desc = "GBP/USD Returns")
getRetNormQuantiles(rnorm(length(Garch$bp)), desc = "Actual Normal")
getRetNormQuantiles(diff(log(Garch$bp)), desc = "log(BP/USD) Returns")
