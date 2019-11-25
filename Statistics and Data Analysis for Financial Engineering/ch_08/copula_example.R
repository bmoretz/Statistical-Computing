library(data.table)
library(ggplot)
library(ggtheme)
library(Ecdat)
library(faraway)
library(fGarch)
library(sn)
library(WVPlots)

theme_set(theme_light())

setwd("D:/Projects/MSDS-RiskAnalytics/datasets")

# Theme Overrides
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             legend.position = "top", legend.title = element_blank())


n <- 10000
data <- data.table( x = runif(n, 0, 1) )

ggplot(data) +
  geom_histogram(aes(x, fill = ..count..))

data$transformed <- qnorm(data$x, 10, 3) 
data$y <- dnorm(data$x, 10, 3)

ggplot(data, aes(x = transformed)) +
  geom_histogram(aes(y = ..density.. , fill = ..count..)) + 
  geom_density(aes( y = ..density.. ), color = "green", lwd = .9) +
  theme(legend.position = "none")

ggplot(data, aes(x, transformed )) +
  geom_line(color = "cornflowerblue", lwd = 1)

ScatterHist(data, "x", "transformed", title = "U -> N")

n <- 10000
df <- 5

data <- data.table( x1 = rt(n, df), x2 = rt(n, df) )

ggplot(data, aes(x = x1)) +
  geom_histogram(aes(y = ..density.. , fill = ..count..)) + 
  geom_density(aes( x = x2, y = ..density.. ), color = "green", lwd = .9) +
  geom_density(aes(x = x2, y = ..density.. ), color = "green", lwd = .9) +
  theme(legend.position = "none")

ggplot(data) +
  geom_density(aes( x = x1, y = ..density.. ), color = "red", lwd = .9) +
  geom_density(aes(x = x2, y = ..density.. ), color = "green", lwd = .9) +
  theme(legend.position = "none")

ScatterHist(data, "x1", "x2", title = "", contour = F)
