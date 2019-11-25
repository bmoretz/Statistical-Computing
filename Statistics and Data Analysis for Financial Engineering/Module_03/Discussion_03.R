library(data.table)
library(ggplot)
library(ggtheme)
library(Ecdat)
library(faraway)
library(fGarch)
library(sn)

theme_set(theme_light())

setwd("D:/Projects/MSDS-RiskAnalytics/datasets")

# Theme Overrides
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             legend.position = "right", legend.title = element_blank())

x <- seq( from = -4, to = 4, by = 0.001)
y_norm <- dnorm(x)
y_df1 <- dt(x, df = 1)
y_df2 <- dt(x, df = 2)
y_df5 <- dt(x, df = 5)
y_df10 <- dt(x, df = 10)
y_df20 <- dt(x, df = 20)


p1 <- ggplot(data = data.table(x = x, y = y_norm)) +
  geom_line(aes(x, y_df1, color = "t-df = 1"), lwd = 1) +
  geom_line(aes(x, y_df2, color = "t-df = 2"), lwd = 1) +
  geom_line(aes(x, y_df5, color = "t-df = 5"), lwd = 1) +
  geom_line(aes(x, y_df10, color = "t-df = 10"), lwd = 1) +
  geom_line(aes(x, y_df20, color = "t-df = 20"), lwd = 1) +
  geom_line(aes(x, y_norm, color = "Normal"), lwd = 1) +
  labs(title = "Normal vs t", x = "", y = "")

p2 <- ggplot(data = data.table(x = x, y = y_norm)) +
  geom_line(aes(x, y_df1, color = "t-df = 1"), lwd = 1) +
  geom_line(aes(x, y_df2, color = "t-df = 2"), lwd = 1) +
  geom_line(aes(x, y_df5, color = "t-df = 5"), lwd = 1) +
  geom_line(aes(x, y_df10, color = "t-df = 10"), lwd = 1) +
  geom_line(aes(x, y_df20, color = "t-df = 20"), lwd = 1) +
  geom_line(aes(x, y_norm, color = "Normal"), lwd = 1) +
  labs(x = "", y = "", title = "Zoomed Left Tail") +
  coord_cartesian(xlim=c(-4, -2), ylim=c(0, .07))

grid.arrange(p1, p2, nrow = 2 )