rm(list = ls())
#---------Fig. 4a--------------

library(vegan)
data <- read.csv("lysogenic_MAG_abundance.csv", row.names=1, header=T)
data <- t(data)
data <- as.data.frame(data)
env <- read.csv("lysogenic_data_env.csv",row.names=1, header=T)
dune_dist <- vegdist(data, method="bray", binary=T) #method="jaccard"
dune_pcoa <- cmdscale(dune_dist, k=3, eig=T)
dune_pcoa_points <- as.data.frame(dune_pcoa$points)
sum_eig <- sum(dune_pcoa$eig)
eig_percent <- round(dune_pcoa$eig/sum_eig*100,1)
colnames(dune_pcoa_points) <- paste0("PCoA", 1:3)
dune_pcoa_result <- cbind(dune_pcoa_points, env)
head(dune_pcoa_result)

library(ggplot2)
ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Country, group=Country)) + #shape=Process
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=4, alpha=0.7
  ) + stat_ellipse(level=0.99) + 
  theme_classic() + scale_shape_manual(values = 15:22)

set.seed(1)

dune.div <- adonis2(data ~ Country, data = env, permutations = 999, method="jaccard", p.adjust.methods="BH")

dune.div

dune_adonis <- paste0("Adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
library(ggplot2)
p1 <- ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=Country)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),title=dune_adonis) +
  geom_point(size=4, alpha=0.7
  ) + scale_color_manual(values = c("China" = "#E3423E","Denmark"= "#2268B1","Argentina"="#38A568", "Singapore"="#5C6AA7", "USA"="#F4A45F")) +
  stat_ellipse(level=0.95) +
  theme_bw() + 
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
        axis.line.x=element_line(linetype="solid",color="black"), #,size=1
        axis.line.y=element_line(linetype="solid",color="black"), #,size=1
        axis.ticks.x=element_line(color="black",linewidth=0.7,lineend = 1.5),
        axis.ticks.y=element_line(color="black",linewidth=0.7,lineend = 1.5),
        axis.title.x=element_text(size=15,family="serif",face="bold",color="black"),
        axis.title.y=element_text(size=15,family="serif",face="bold",color="black"),
        plot.title = element_text(size=14,family="serif",face="bold",color="black", hjust = 0.04, vjust = 0,margin=margin(t=0,b=-14)),
        axis.text.x = element_text(size=14,family="serif",face="bold",color="black"),
        axis.text.y = element_text(size=14,family="serif",face="bold",color="black"),
        legend.title=(element_text(size=14,family="serif",face="bold",color="black")),
        legend.background =element_blank(),
        legend.text = (element_text(size=14,family="serif",face="bold",color="black")),
        legend.position = c(0.8,0.8),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

p1

#--------Fig. 4b----------
rm(list = ls())

library(vegan)
votu<-read.csv('prophage_abundance.csv',header = T, row.names = 1)
votu<-decostand(t(votu),method='hellinger')
distv<-vegdist(votu,method = 'bray')
distv_num<-1-as.numeric(distv)

library(SoDA)
loca<-read.csv('site.lat_lng.csv',header = T, row.names = 1)
xy_loca<-geoXY(loca$lat,loca$lng)
rownames(xy_loca)<-rownames(loca)
dist_loca<-vegdist(xy_loca,method = 'euclidean')
dist_loca_log<-log10(dist_loca)
dist_loca_num<-as.numeric(dist_loca_log)

data_new<-data.frame(distv_num,dist_loca_num)
data_new <- subset(data_new, dist_loca_num > 0)
data_new$group<-ifelse(data_new$dist_loca_num < 4, "Local","Regional")
data_new$group<-as.factor(data_new$group)

data_local<-data_new[which(data_new$group == "Local"),]
data_regional<-data_new[which(data_new$group == "Regional"),]

library(ggplot2)
p1 <- ggplot(data_new)+geom_point(size=1,color="grey",aes(x = dist_loca_num, y = distv_num))+
  labs(x="Log10 geographic distance (m)",y="prophage community taxonomic similarity")+
  geom_smooth(method="lm",fill= NA,level=0,aes(dist_loca_num,distv_num,color=group))+
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
        panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15,family="serif",face="bold",color="black"),
        axis.title.y = element_text(size=15,family="serif",face="bold",color="black"),
        axis.text.x = element_text(hjust =0.5,size=14,family="serif",face="bold",color="black"),
        axis.text.y=element_text(size=14,family="serif",face="bold",color="black"),
        legend.text = element_text(size=14))+expand_limits(x=0,y=0)+
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1))+
  scale_x_continuous(breaks = seq(0,7,1))+
  scale_shape_manual(values=c(16,16))+
  scale_color_manual(values=c("#E3423E", "#2268B1", "#38A568"))+
  geom_smooth(method="lm",fill= NA,level=0,aes(dist_loca_num,distv_num,color="#E3423E"))+
  theme(legend.position = c(0.85,0.9), legend.background = element_blank())#+coord_cartesian(ylim=c(0,1),xlim = c(0,10))

p1

#-------Fig. 4c-------------
rm(list = ls())

library(UpSetR)

otu <- read.csv("lysogenic_MAG_Country_sum.csv", header = TRUE, row.names = 1)
otu1 <- otu
otu1[otu1 > 0] <- 1

upset(otu1, nset = 6, nintersects = 100, mainbar.y.label = "Numbers of lysogenic MAG", 
      mainbar.y.max = 270, sets.x.label = "Numbers of lysogenic MAG",
      order.by = c('freq'), decreasing = c(TRUE, TRUE))