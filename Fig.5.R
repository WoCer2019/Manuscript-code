#------------Fig.5c--------------

mydata <- read.csv("induction.csv", header = T)
cbPalette <- c("#2268B1","#38A568", "#5C6AA7", "#F4A45F","#E3423E","#E69F00", "#CC79A7", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#999999","#0072B2","#D55E00")

library(ggplot2)
library(gghalves)

compaired <- list(c("CK", "MMC"))
(p <- ggplot(mydata,aes(x=group, y=Num, color=group, fill=group))+
    scale_fill_manual(values = cbPalette)+
    scale_colour_manual(values = cbPalette)+
    geom_half_violin(position=position_nudge(x=0.15,y=0),
                     side='R',adjust=1.5,trim=F,alpha=0.4,colour = NA)+
    geom_half_point(position=position_nudge(x=-0.43,y=0),size =3, shape =19,range_scale = 0.5,alpha=0.5)+
    geom_boxplot(alpha=0.4,width=0.2,
                 position=position_dodge(width=0.8),
                 size=0.75,outlier.colour = NA)+
    theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), #linewidth=1
          axis.ticks.x=element_line(color="black",size=1,lineend = 10),
          axis.ticks.y=element_line(color="black",size=1,lineend = 10),
          axis.title.x=element_text(size=15,family="serif",face="bold",color="black"),
          axis.title.y=element_text(size=15,family="serif",face="bold",color="black"),
          plot.title = element_text(size=14,family="serif",face="bold",color="black", hjust = 0.04, vjust = 0,margin=margin(t=0,b=-14)),
          axis.text.x = element_text(size=14,family="serif",face="bold",color="black"), #, angle = 0, hjust = 1, vjust = 1
          axis.text.y = element_text(size=14,family="serif",face="bold",color="black"),
          legend.title=(element_text(size=14,family="serif",face="bold",color="black")),
          legend.background =element_blank(),
          legend.text = (element_text(size=14,family="serif",face="bold",color="black")),
          legend.position = c(0.88,0.1),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()
    ) + 
    geom_signif(comparisons = compaired,
                step_increase = 0.1,
                map_signif_level = F,
                test = wilcox.test, 
                textsize = 4.5, family="serif", color="#000000",
                size = 0.5) +
    ylab("Virus-like particels /mL") +
    xlab(NULL)
)

#------------Fig.5d

#------------Fig.5d----------
library(ggplot2)
library(dplyr)
library(ggrepel)

rm(list = ls())

df <- read.csv("abundance_decreasing_MAG.csv", sep = ",", stringsAsFactors = FALSE)
df$NegativeLogP=-log10(df$p_value)

df$MAG <- factor(df$MAG, levels = df$MAG[order(df$group)])

getBreak<-function(x,y){
  freq<-as.vector(table(y)) 
  half_freq<-freq%/%2
  for (i in seq(2,length(freq))){
    new_num<-freq[i]+freq[i-1]
    freq[i]=new_num}
  pos<-freq-half_freq
  break_point<-as.vector(x[pos])
  return(break_point)
}

gtext<-df %>% filter(MAG %in% c("AS_CK_bin.28", "AS_MMC_bin.16","AS_MMC_bin.52","AS_MMC_bin.56","AS_MMC_bin.31","AS_MMC_bin.188"))


threshold<--log10(0.05)

mycol<-c("#264654","#8AB07C","#299D91","#86C03D",
         "#E8C56A","#F3A263","#E56F51","#AF4592",
         "#7AC5C1","#906024","#6E1F86","#545556")


(p<-ggplot(df, aes(x=MAG, y=NegativeLogP, color=group, size=Decreasing.percentage.of.abundance, shape=enrich)) +
    annotate("rect",xmin = 1, xmax = 16,
             ymin = 0, ymax = 5,
             fill="#F4F4F3")+
    annotate("rect",xmin = 25, xmax = 51,
             ymin = 0, ymax = 5,
             fill="#F4F4F3")+
    annotate("rect",xmin = 70, xmax = 87,
             ymin = 0, ymax = 5,
             fill="#F4F4F3")+
    annotate("rect",xmin = 96, xmax = 117,
             ymin = 0, ymax = 5,
             fill="#F4F4F3")+
    geom_point(alpha=.7,key_glyph="point")+
    scale_shape_manual(values=c(1,16))+
    guides(color="none",shape="none")+
    scale_x_discrete(breaks=getBreak(df$MAG,df$group),
                     labels=c("Acidobacteriota","Actinobacteriota","Bacteroidota","Chloroflexota","Others","Planctomycetota","Planctomycetes"),
                     expand = expansion(mult = 0.03))+
    geom_text_repel(aes(MAG,NegativeLogP,label=c("Skermania","unclassified Xanthobacteraceae","Azonexus",
                                                 "Ca. Accumulibacter","unclassified Burkholderiaceae","Bradyrhizobium")),
                    gtext,colour='black', size =4,box.padding = 3, 
                    point.padding = 0.8, 
                    segment.color = "black",
                    show.legend = F)+
    geom_hline(yintercept=threshold, linetype=2, color="grey")+
    labs(x="", y="-log10(P)")+ 
    scale_color_manual(values = mycol)+ 
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), #linewidth=1
          axis.ticks.x=element_line(color="black",size=1,lineend = 10),
          axis.ticks.y=element_line(color="black",size=1,lineend = 10),
          axis.title.x=element_text(size=15,family="serif",face="bold",color="black"),
          axis.title.y=element_text(size=15,family="serif",face="bold",color="black"),
          plot.title = element_text(size=14,family="serif",face="bold",color="black", hjust = 0.04, vjust = 0,margin=margin(t=0,b=-14)),
          axis.text.x = element_text(size=14,family="serif",face="bold",color="black", angle = -45, hjust = 0.01, vjust = 1),
          axis.text.y = element_text(size=14,family="serif",face="bold",color="black"),
          legend.title=(element_text(size=14,family="serif",face="bold",color="black")),
          legend.background =element_blank(),
          legend.text = (element_text(size=14,family="serif",face="bold",color="black")),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_rect(fill = 'white'))
)

#----------Fig.5ef---------
rm(list = ls())
df <- read.csv("active_lysogenic_function_MAG_country_abundance_sum.csv", 
               header = T)
df[df == 0] <- NA
library(reshape)
dat<-reshape::melt(df,id='X')
dat$X <- factor(dat$X, levels = unique(dat$X))
library(ggplot2)
(p <- ggplot(dat, aes(variable, X)) +
    geom_point(aes(size = value, color=variable), shape = 16) +
    scale_color_manual(values = c("China" = "#E3423E","Denmark"= "#2268B1",
                                  "Argentina"="#38A568", "Singapore"="#5C6AA7", 
                                  "USA"="#F4A45F"))+
    scale_size_continuous(range = c(0.1, 5)) +
    labs(x = "", y = "", title = "") +
    theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          axis.ticks.x=element_line(color="black",size=1,lineend = 10),
          axis.ticks.y=element_line(color="black",size=1,lineend = 10),
          axis.title.x=element_text(size=15,family="serif",face="bold",color="black"),
          axis.title.y=element_text(size=15,family="serif",face="bold",color="black"),
          plot.title = element_text(size=14,family="serif",face="bold",color="black", hjust = 0.04, vjust = 0,margin=margin(t=0,b=-14)),
          axis.text.x = element_text(size=14,family="serif",face="bold",color="black", angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size=14,family="serif",face="bold",color="black"),
          legend.title=(element_text(size=14,family="serif",face="bold",color="black")),
          legend.background =element_blank(),
          legend.text = (element_text(size=14,family="serif",face="bold",color="black")),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()
    )
)