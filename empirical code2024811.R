######################################selecting dataset
library(rtry);library(patchwork)
daf=rtry_import("C:\\Users\\Lenovo\\Desktop\\35431.txt")

DF=rtry_exclude(daf,
             (TraitID %in% NA),
             baseOn = ObsDataID  )
DF2=rtry_exclude(DF,
                (DataID %in% c(8143,8144,8273,910,1361,1328,1329,1330,2034,472,573,887,899,378)),
                baseOn = ObsDataID)

#The VISTA Plant Trait Database
DF4=rtry_select_row(daf,
                    (DatasetID %in% 45),
                    getAncillary = TRUE)

traits=rtry_select_anc(DF4,114,62,59,60,80,319,92,490,
                      491,492,493,271,454,494,495,496,
                      497,498,499,500,501,502,61)

data_selected <- rtry_select_col(DF4,
                                 DatasetID,	Dataset,	SpeciesName,	AccSpeciesID,	AccSpeciesName,	ObservationID,
                                 ObsDataID,	TraitID,	TraitName)

data_selected2=data_selected[!duplicated(data_selected$ObservationID),]
all=rtry_join_left(data_selected2,traits,baseOn = ObservationID)

#write.csv(DF2,"flowering.csv")
#write.csv(DF4,"grazingeffect.csv")
#write.csv(all,"gr2.csv")
######################################
library(ggplot2)
library(rtry)
df=read.csv("flowering.csv")
mm=sort(table(df$SpeciesName),decreasing = TRUE)#
df2=rtry_select_row(df,
                    (SpeciesName %in% c("Brassica tournefortii")),
                    getAncillary = TRUE)

df2$OrigValueStr <- as.factor(df2$OrigValueStr)
freq=table(df2$OrigValueStr)
df3=data.frame(freq=freq,NUM=1:length(freq))
figempA=ggplot(df2, aes(x=OrigValueStr)) + geom_bar(fill="steelblue")+  
  geom_smooth(data = df3,aes(x=NUM,y=freq.Freq),method="loess",span=0.35,size=2,se=FALSE)+
  labs(title="Multiple stable distribution",colour="",x="Flowering time (day)",y=expression(paste("Frequency ")))+
  annotate("text",x=1.5, y=230, label="A",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=15,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.80,hjust=0.8,angle = 45),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_y_continuous(breaks=seq(0,240,by=60),limits=c(0,240))



#########effect of disturbance (e.g. grazing) on trait evolution
df <- data.frame( category = c("significant", "not significant"),
  value = c(82, 18) )
figempB=ggplot(df, aes(x="", y=value,fill=category)) +  geom_bar(stat="identity")+coord_polar("y", start=0)+
  labs(title="Disturbance induced \n evolution",colour="",x="",y=expression(paste("")))+
  annotate("text",x=1.2, y=90, label="No significant \n (18%)",colour="black",angle = 40,size=6 ,family="serif")+
  annotate("text",x=1, y=50, label="Significant (82%)",colour="black",angle = 0,size=6,family="serif")+
  theme(axis.text=element_text(size=15,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.80,hjust=0.8,angle = 0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = "none", 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )
#####################
df=read.csv("gr2.csv")# analyses about the relationship of biomass vs. onset of flowering; 
df2=rtry_select_row(df,
                    (SpeciesName %in% c("Plantago lanceolata","Dactylis glomerata","Trifolium pratense",
                                        "Anthoxanthum odoratum","Ranunculus acris")),
                    getAncillary = TRUE)
figempC=ggplot(data = df2,aes(x=Onsetflower,y=AGBmax,group=SpeciesName,color=SpeciesName))+geom_point(size=4)+geom_line(size=2)+
  labs(title="Intraspecies level",colour="",x="Flowering time (Julian day)",y=expression(paste("Vegetation biomass ", (g/m^2))))+
  annotate("text",x=60, y=720, label="C",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.82,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(50,250,by=50),limits=c(50,250))+
  scale_y_continuous(breaks=seq(50,750,by=175),limits=c(50,750))
##############insterspecies level
df2=rtry_select_row(df,
                    (Genus %in% c("Trifolium")),
                    getAncillary = TRUE)
meanAGBmax=aggregate(df2$AGBmax,by=list(type=df2$Onsetflower),mean)
bioup=Onflowerup=biodwn=Onflowerdwn=NA;j=k=1
for (i in 1:length(meanAGBmax$type)){
  if(meanAGBmax$x[i]>=250){bioup[j]=meanAGBmax$x[i];Onflowerup[j]=meanAGBmax$type[i];j=j+1}
  else {biodwn[k]=meanAGBmax$x[i];Onflowerdwn[k]=meanAGBmax$type[i];k=k+1}
}
dfup=data.frame(bio=bioup,Onflower=Onflowerup)
dfdwn=data.frame(bio=biodwn,Onflower=Onflowerdwn)
figempD=ggplot()+
  geom_point(data = dfup,aes(x=Onflower,y=bio),size=4,colour="orange")+geom_smooth(data = dfup,aes(x=Onflower,y=bio),method="loess",span=1,size=2)+
  geom_point(data = dfdwn,aes(x=Onflower,y=bio),size=4,colour="orange")+geom_smooth(data = dfdwn,aes(x=Onflower,y=bio),method="loess",span=1,size=2)+
  labs(title="Interspecies level (Trifolium)",colour="",x="Flowering time (Julian day)",
       y=expression(paste("Vegetation biomass ", (g/m^2))))+
  annotate("text",x=55, y=720, label="D",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(50,250,by=50),limits=c(50,250))+
  scale_y_continuous(breaks=seq(50,750,by=175),limits=c(50,750))

#####################appendix figure
df2=rtry_select_row(df,
                    (SpeciesName %in% c("Plantago lanceolata")),
                    getAncillary = TRUE)

df2$Onsetflower <- as.factor(df2$Onsetflower)
freq=table(df2$Onsetflower)
df3=data.frame(freq=freq,NUM=1:length(freq))
figempS1A=ggplot(df2, aes(x=Onsetflower)) + geom_bar(fill="steelblue")+  
  geom_smooth(data = df3,aes(x=NUM,y=freq.Freq),method="loess",span=0.35,size=2,se=FALSE)+
  labs(title="Plantago lanceolata",colour="",x="Flowering time (Julian day)",y=expression(paste("Frequency ")))+
  annotate("text",x=1.5, y=3.6, label="A",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=15,family="serif",colour = "black"),
        axis.text.x = element_text(angle = 0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_y_continuous(breaks=seq(0,4,by=1),limits=c(0,4))


df2=rtry_select_row(df,
                    (SpeciesName %in% c("Dactylis glomerata")),
                    getAncillary = TRUE)

df2$Onsetflower <- as.factor(df2$Onsetflower)
freq=table(df2$Onsetflower)
df3=data.frame(freq=freq,NUM=1:length(freq))
figempS1B=ggplot(df2, aes(x=Onsetflower)) + geom_bar(fill="steelblue")+  
  geom_smooth(data = df3,aes(x=NUM,y=freq.Freq),method="loess",span=0.35,size=2,se=FALSE)+
  labs(title="Dactylis glomerata",colour="",x="Flowering time (Julian day)",y=expression(paste("Frequency ")))+
  annotate("text",x=1.5, y=3.6, label="B",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=15,family="serif",colour = "black"),
        axis.text.x = element_text(),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_y_continuous(breaks=seq(0,4,by=1),limits=c(0,4))


df2=rtry_select_row(df,
                    (SpeciesName %in% c("Trifolium pratense")),
                    getAncillary = TRUE)

df2$Onsetflower <- as.factor(df2$Onsetflower)
freq=table(df2$Onsetflower)
df3=data.frame(freq=freq,NUM=1:length(freq))
figempS1C=ggplot(df2, aes(x=Onsetflower)) + geom_bar(fill="steelblue")+  
  geom_smooth(data = df3,aes(x=NUM,y=freq.Freq),method="loess",span=0.35,size=2,se=FALSE)+
  labs(title="Trifolium pratense",colour="",x="Flowering time (Julian day)",y=expression(paste("Frequency ")))+
  annotate("text",x=1.5, y=3.6, label="C",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=15,family="serif",colour = "black"),
        axis.text.x = element_text(),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_y_continuous(breaks=seq(0,4,by=1),limits=c(0,4))


df2=rtry_select_row(df,
                    (SpeciesName %in% c("Anthoxanthum odoratum")),
                    getAncillary = TRUE)

df2$Onsetflower <- as.factor(df2$Onsetflower)
freq=table(df2$Onsetflower)
df3=data.frame(freq=freq,NUM=1:length(freq))
figempS1D=ggplot(df2, aes(x=Onsetflower)) + geom_bar(fill="steelblue")+  
  geom_smooth(data = df3,aes(x=NUM,y=freq.Freq),method="loess",span=0.35,size=2,se=FALSE)+
  labs(title="Anthoxanthum odoratum",colour="",x="Flowering time (Julian day)",y=expression(paste("Frequency ")))+
  annotate("text",x=1.5, y=3.6, label="D",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=15,family="serif",colour = "black"),
        axis.text.x = element_text(),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_y_continuous(breaks=seq(0,4,by=1),limits=c(0,4))

df2=rtry_select_row(df,
                    (SpeciesName %in% c("Ranunculus acris")),
                    getAncillary = TRUE)

df2$Onsetflower <- as.factor(df2$Onsetflower)
freq=table(df2$Onsetflower)
df3=data.frame(freq=freq,NUM=1:length(freq))
figempS1E=ggplot(df2, aes(x=Onsetflower)) + geom_bar(fill="steelblue")+  
  geom_smooth(data = df3,aes(x=NUM,y=freq.Freq),method="loess",span=0.35,size=2,se=FALSE)+
  labs(title="Ranunculus acris",colour="",x="Flowering time (Julian day)",y=expression(paste("Frequency ")))+
  annotate("text",x=1.2, y=3.6, label="E",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=15,family="serif",colour = "black"),
        axis.text.x = element_text(),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_y_continuous(breaks=seq(0,4,by=1),limits=c(0,4))


##########interspecies
df2=rtry_select_row(df,
                    (Genus %in% c("Carex")),
                    getAncillary = TRUE)

df2$Onsetflower <- as.factor(df2$Onsetflower)
freq=table(df2$Onsetflower)
df3=data.frame(freq=freq,NUM=1:length(freq))
figempS2A=ggplot(df2, aes(x=Onsetflower)) + geom_bar(fill="steelblue")+  
  geom_smooth(data = df3,aes(x=NUM,y=freq.Freq),method="loess",span=0.35,size=2,se=FALSE)+
  labs(title="Interspecies level (Carex)",colour="",x="Flowering time (Julian day)",y=expression(paste("Frequency ")))+
  annotate("text",x=1.5, y=14, label="A",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=15,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.80,hjust=0.8,angle = 45),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_y_continuous(breaks=seq(0,15,by=5),limits=c(0,15))

#########
meanAGBmax=aggregate(df2$AGBmax,by=list(type=df2$Onsetflower),mean)
meanAGBmax=meanAGBmax[complete.cases(meanAGBmax),]#
bioup=Onflowerup=biodwn=Onflowerdwn=NA;j=k=1
for (i in 1:length(meanAGBmax$type)){
  if(meanAGBmax$x[i]>=200){bioup[j]=meanAGBmax$x[i];Onflowerup[j]=meanAGBmax$type[i];j=j+1}
  else {biodwn[k]=meanAGBmax$x[i];Onflowerdwn[k]=meanAGBmax$type[i];k=k+1}
}
dfup=data.frame(bio=bioup,Onflower=Onflowerup)
dfdwn=data.frame(bio=biodwn,Onflower=Onflowerdwn)

figempS2B=ggplot()+
  geom_point(data = dfup,aes(x=Onflower,y=bio),size=4,colour="orange")+geom_smooth(data = dfup,aes(x=Onflower,y=bio),method="loess",span=1,size=2)+
  geom_point(data = dfdwn,aes(x=Onflower,y=bio),size=4,colour="orange")+geom_smooth(data = dfdwn,aes(x=Onflower,y=bio),method="loess",span=1,size=2)+
  labs(title="Interspecies level (Carex)",colour="",x="Flowering time (Julian day)",
       y=expression(paste("Vegetation biomass ", (g/m^2))))+
  annotate("text",x=0, y=600, label="B",colour="black",angle = 0,size=10 ,family="serif")+
  theme(axis.text=element_text(size=30,family="serif",colour = "black"),
        axis.text.x = element_text(vjust=0.0),
        axis.title=element_text(size=30,family="serif", angle=0, face="plain"),
        axis.title.x = element_text(vjust=0.0),
        axis.line=element_line(size = 0.1, colour = "black", linetype=1),
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(size = 0.5, color="black"),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(size = 2.0,colour = "black",fill = "NA", linetype=1),
        panel.background = element_rect(fill = "white", color = NA),
        line = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.75,0.8), 
        legend.background = element_rect(fill = "NA", color = NA),
        legend.text =element_text(size = 16,color='black'),
        legend.title = element_text(size=16, face = "bold"),
        legend.direction="vertical", #horizontal; vertical
        legend.key.size = unit(2,'line'), # space between legend text
        plot.title = element_text(size=30, hjust=0.5, color = "black", face = "bold",
                                  margin = margin(b = -0.0, t = 0.4, l = 0, unit = "cm")),
        plot.margin = margin(t=0.2, r=0.3, b=0.2, l=0.12, "cm") )+
  scale_x_continuous(breaks=seq(0,40,by=10),limits=c(0,40))+
  scale_y_continuous(breaks=seq(0,600,by=200),limits=c(0,600))
####################################
#ggsave("figempirical.pdf",(figempA|figempB)/(figempC|figempD),width = 50, height = 35, units = "cm", dpi = 300) 
#ggsave("figSintrafreq.pdf",(figempS1A|figempS1B)/(figempS1C|figempS1D)/figempS1E,width = 40, height = 45, units = "cm", dpi = 300) 
#ggsave("figSinter.pdf",figempS2A/figempS2B,width = 30, height = 30, units = "cm", dpi = 300) 

####################################




