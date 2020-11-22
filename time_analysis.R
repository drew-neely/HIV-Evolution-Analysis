
library(ggplot2)
library(tidyverse)


data <- read.table("time_distance_matrix.csv", header=T, dec=".", sep=",")
data$patient <- as.factor(data$patient)

ggplot(data, aes(x=days, y=dist, color=patient)) +
  geom_point() +
  geom_smooth(method='loess', se=FALSE) +
  coord_cartesian(xlim=c(10,265), ylim=c(0.0002,0.007)) + 
  ggtitle("Genetic Distance in HIV Genome over Time", subtitle = "Measured from first sample") +
  xlab("Days since seroconversion") + 
  ylab("Genetic distance from first sample (GGDC formula 2 distance)") +
  theme(plot.title = element_text(hjust = 0.5,size = 14), plot.subtitle = element_text(hjust = 0.5,size=10), 
        axis.text.y = element_text(angle=50, hjust=0.5), 
        axis.title = element_text(size = 13))


