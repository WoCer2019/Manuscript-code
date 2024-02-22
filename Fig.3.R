rm(list = ls())
library(ggplot2)
library(gg.gap)
windowsFonts(A=windowsFont("Times New Roman"),
             B=windowsFont("Arial"))
theme(text=element_text(family="A",face="bold",size=20))


data <- read.csv("viral_contigs_length_distribution.csv",header = T)
data_part <- read.csv("viral_contigs_length_distribution.csv",header = T)


p <- ggplot(data=data)+geom_histogram(mapping = aes(x = length), binwidth = 4, fill="#2268B1", alpha=0.6)+
  xlab("The numbers of prophage in each lysogenic MAG") + ylab("Count") +
  theme_bw() + 
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
        axis.title.x=element_text(size=15,family="serif",face="bold",color="black"),
        axis.title.y=element_text(size=15,family="serif",face="bold",color="black"),
        axis.text.x = element_text(size=14,family="serif",face="bold",color="black"),
        axis.text.y = element_text(size=14,family="serif",face="bold",color="black"),
        legend.title=(element_text(size=14,family="serif",face="bold",color="black")),
        legend.background =element_blank(),
        legend.text = (element_text(size=14,family="serif",face="bold",color="black")),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
p

p_part <- ggplot(data=data_part)+geom_histogram(mapping = aes(x = length), binwidth = 1, fill="#2268B1", alpha=0.6)+
  xlab("The numbers of prophage in each lysogenic MAG") + ylab("Count") + theme(axis.line = element_line(colour = "black")) +
  theme_classic() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(#panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
    axis.title.x=element_text(size=15,family="serif",face="bold",color="black"),
    axis.title.y=element_text(size=15,family="serif",face="bold",color="black"),
    axis.text.x = element_text(size=14,family="serif",face="bold",color="black"),
    axis.text.y = element_text(size=14,family="serif",face="bold",color="black"),
    legend.title=(element_text(size=14,family="serif",face="bold",color="black")),
    legend.background =element_blank(),
    legend.text = (element_text(size=14,family="serif",face="bold",color="black")),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  )

p_part