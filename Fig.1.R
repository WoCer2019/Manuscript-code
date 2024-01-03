rm(list = ls())
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(dplyr)
library(ggtrendline)
data <- read.csv("fit.csv",header = T)
x <- data$Completeness
y <- data$numbers


(p <- ggtrendline(x,y,"exp3P",linecolor = "black",linetype = "solid",linewidth = 1,CI.level = 0.95,
                  CI.fill = "grey60",CI.alpha = 0.3,CI.color = "black",CI.lty = 0,CI.lwd = 0.5,    #CI.lty	line type of confidence interval. CI.lwd	line width of confidence interval.
                  eq.x = NULL, eq.y = NULL, rrp.x =NULL, rrp.y =NULL, #equation position, the position for R square and P value
                  Rname = 1,Pname = 0, eDigit = 2, eSize = 3) + 
  geom_point(aes(x=x, y=y),color = "#33A02C", size = 3) +
  labs(x="MAG completeness (%)",y = "Percentage of lysogenic MAGs in total MAGs (%)", title= "") +
  theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), #linewidth=1
        # axis.line.x=element_line(linetype="solid",color="black"), # size=1
        # axis.line.y=element_line(linetype="solid",color="black"),
        axis.ticks.x=element_line(color="black",size=1,lineend = 10),
        axis.ticks.y=element_line(color="black",size=1,lineend = 10),
        axis.title.x=element_text(size=15,family="serif",face="bold",color="black"),
        axis.title.y=element_text(size=15,family="serif",face="bold",color="black"),
        plot.title = element_text(size=14,family="serif",face="bold",color="black", hjust = 0.04, vjust = 0,margin=margin(t=0,b=-14)),
        axis.text.x = element_text(size=14,family="serif",face="bold",color="black"), #, angle = 45, hjust = 1, vjust = 1
        axis.text.y = element_text(size=14,family="serif",face="bold",color="black"),
        legend.title=(element_text(size=14,family="serif",face="bold",color="black")),
        legend.background =element_blank(),
        legend.text = (element_text(size=14,family="serif",face="bold",color="black")),
        legend.position = c(0.88,0.1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
  ))