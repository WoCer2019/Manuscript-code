rm(list = ls())
#---------Fig. 4--------------
library(UpSetR)

otu <- read.csv("lysogenic_MAG_Country_sum.csv", header = TRUE, row.names = 1)
otu1 <- otu
otu1[otu1 > 0] <- 1

upset(otu1, nset = 6, nintersects = 100, mainbar.y.label = "Numbers of lysogenic MAG", 
      mainbar.y.max = 270, sets.x.label = "Numbers of lysogenic MAG", #mainbar.y.max = 620
      order.by = c('freq'), decreasing = c(TRUE, TRUE))
