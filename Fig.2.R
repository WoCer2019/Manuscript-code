#----------Fig.2b-------------------------

rm(list = ls())
mydata<-read.csv("Phylum_lysogen_percentage_based_numbers.csv",sep=",",na.strings="NA",stringsAsFactors=FALSE)

mydata$Phylum<- factor(mydata$Phylum, levels = c("Nitrospirota (28 MAGs)","Planctomycetota (75 MAGs)","Myxococcota (80 MAGs)","Acidobacteriota (103 MAGs)",
                                                 "Proteobacteria (427 MAGs)","Bacteroidota (393 MAGs)","Chloroflexota (136 MAGs)",
                                                 "Verrucomicrobiota (34 MAGs)","Gemmatimonadota (30 MAGs)","Actinobacteriota (125 MAGs)","Patescibacteria (202 MAGs)"))

(ggplot(mydata, aes(Phylum, Percentage)) +
    geom_segment(aes(x=Phylum, 
                     xend=Phylum,
                     y=55, 
                     yend=Percentage),linewidth=0.7)+ 
    geom_point(shape=21,size=3,colour="black",fill="#2268B1")+
    geom_hline(yintercept = 55, linetype = "dashed", color = "#3B9D66", linewidth=1)+
    theme_bw() + 
    labs(y = "Percentage of lysogenic MAGs (%)", x = NULL) +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"), #linewidth=1
          axis.ticks.x=element_line(color="black",linewidth=1,lineend = 10),
          axis.ticks.y=element_line(color="black",linewidth=1,lineend = 10),
          axis.title.x=element_text(size=14,family="serif",face="bold",color="black"),
          axis.title.y=element_text(size=14,family="serif",face="bold",color="black"),
          plot.title = element_text(size=14,family="serif",face="bold",color="black", hjust = 0.04, vjust = 0,margin=margin(t=0,b=-14)),
          axis.text.x = element_text(size=13,family="serif",face="bold",color="black", angle = -30, hjust = 0, vjust = 1),
          axis.text.y = element_text(size=13,family="serif",face="bold",color="black"),
          legend.title=(element_text(size=14,family="serif",face="bold",color="black")),
          legend.background =element_blank(),
          legend.text = (element_text(size=14,family="serif",face="bold",color="black")),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()
    ) -> p3)

#----------Fig.2c---------------
mydata<-read.csv("Phylum_lysogen_percentage_based_abundance.csv",sep=",",na.strings="NA",stringsAsFactors=FALSE)

mydata$Phylum<- factor(mydata$Phylum, levels = c("Gemmatimonadota (30 MAGs)","Nitrospirota (28 MAGs)","Chloroflexota (136 MAGs)",
                                                 "Myxococcota (80 MAGs)","Verrucomicrobiota (34 MAGs)","Bacteroidota (393 MAGs)",
                                                 "Proteobacteria (427 MAGs)","Acidobacteriota (103 MAGs)",
                                                 "Actinobacteriota (125 MAGs)","Planctomycetota (75 MAGs)","Patescibacteria (202 MAGs)"))

(ggplot(mydata, aes(Phylum, Percentage)) +
    geom_segment(aes(x=Phylum, 
                     xend=Phylum,
                     y=59, 
                     yend=Percentage),linewidth=0.7)+
    geom_point(shape=21,size=3,colour="black",fill="#2268B1")+
    geom_hline(yintercept = 59, linetype = "dashed", color = "#3B9D66", linewidth=1)+
    theme_bw() +
    labs(y = "Relative abundance of lysogenic MAGs (%)", x = NULL) +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"), #linewidth=1
          axis.ticks.x=element_line(color="black",linewidth=1,lineend = 10),
          axis.ticks.y=element_line(color="black",linewidth=1,lineend = 10),
          axis.title.x=element_text(size=13,family="serif",face="bold",color="black"),
          axis.title.y=element_text(size=13,family="serif",face="bold",color="black"),
          plot.title = element_text(size=13,family="serif",face="bold",color="black", hjust = 0.04, vjust = 0,margin=margin(t=0,b=-14)),
          axis.text.x = element_text(size=12,family="serif",face="bold",color="black", angle = -90, hjust = 0, vjust = 1),
          axis.text.y = element_text(size=12,family="serif",face="bold",color="black"),
          legend.title=(element_text(size=13,family="serif",face="bold",color="black")),
          legend.background =element_blank(),
          legend.text = (element_text(size=13,family="serif",face="bold",color="black")),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()
    ) -> p3 )