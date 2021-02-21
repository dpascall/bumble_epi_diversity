library(MCMCglmm)
library(ggplot2)
library(prevalence)
library(SPIn)
rm(list = ls())

data<-read.csv("~/Dropbox/P1datafinal.csv")[c(1:759),]
data<-data[,c("Species","Location","RLV","LMV","MV1","MV2","SBPV","ABPV")]
colnames(data)<-c("Species","Location","River Luinaeg virus","Loch Morlich virus","Mayfield virus 1","Mayfield virus 2","SBPV","ABPV")
reshapeddata<-reshape(data, direction="long", varying=list(names(data)[-c(1:2)]), times=names(data)[-c(1:2)])[,c(1:4)]
colnames(reshapeddata)<-c("Host","Location","Virus","Status")
reshapeddata$Host<-as.character(reshapeddata$Host)
for (i in 1:length(unique(reshapeddata$Host))) {
  reshapeddata[which(reshapeddata$Host==unique(reshapeddata$Host)[i]),"Host"]<-paste("B.",unique(reshapeddata$Host)[i],sep=" ")
}

reshapeddata$Host<-as.factor(reshapeddata$Host)

library(prevalence)
library(SPIn)

mat<-matrix(NA,nrow=length(unique(reshapeddata$Species))*length(unique(reshapeddata$Location))*length(unique(reshapeddata$Virus)),ncol=8)
colnames(mat)<-c("Virus","Location","Host","LL","L","M","H","HH")
mat<-as.data.frame(mat)

k<-1
for (i in 1:length(unique(reshapeddata$Virus))) {
  for (j in 1:length(unique(reshapeddata$Location))) {
    for (u in 1:length(unique(reshapeddata$Host))) {
      temp<-subset(reshapeddata,Virus==unique(reshapeddata$Virus)[i]&Host==unique(reshapeddata$Host)[u]&Location==unique(reshapeddata$Location)[j])
      temp<-temp[!(temp$Status%in%NA),]
      mat[k,1]<-as.character(unique(reshapeddata$Virus)[i])
      mat[k,2]<-as.character(unique(reshapeddata$Location)[j])
      mat[k,3]<-as.character(unique(reshapeddata$Host)[u])
      if (nrow(temp)>0){
        message(nrow(temp))
        sims<-rbeta(20000,(1+sum(temp$Status)),(1+(nrow(temp)-sum(temp$Status))))
        mat[k,6]<-posterior.mode(as.mcmc(sims))
        mat[k,c(4,8)]<-SPIn(as.mcmc(sims),conf=0.9,lb=0,ub=1)$spin
        mat[k,c(5,7)]<-SPIn(as.mcmc(sims),conf=0.5,lb=0,ub=1)$spin
        print(k)
      }
      k<-k+1
    }
  }
}  

mat

library(ggplot2)

colours <- c("#68023F", "#008169", "#EF0096", "#00DCB5", "#FFCFE2", "#003C86", "#9400E6", "#009FFA",
             "#FF71FD", "#7CFFFA", "#6A0213", "#008607", "#00E307")

coeffplot <- ggplot(mat)
coeffplot <- coeffplot + facet_grid(Location~Virus)
coeffplot <- coeffplot + geom_point(aes(x = Host, y = M, colour = Host), position = position_dodge(width = 1/2), size=2)
coeffplot <- coeffplot + geom_linerange(aes(x = Host, ymin = L,ymax = H, colour = Host),lwd = 2, position = position_dodge(width = 1/2))
coeffplot <- coeffplot + geom_linerange(aes(x = Host, ymin = LL,ymax = HH, colour = Host),lwd = 1, position = position_dodge(width = 1/2))
coeffplot <- coeffplot + theme_bw() + scale_color_manual(values = colours) + theme(legend.position = "none")
coeffplot <- coeffplot + ylab("Prevalence") + theme(axis.text.x = element_text(angle = 45, hjust = 1,face = "italic", size = 5.5))
print(coeffplot)  # The trick to these is position_dodge()

tiff("~/Desktop/HostByLoc.tif", width = 189, height = 267.3, res = 300, units = "mm")
coeffplot
dev.off()

##stacked bar chart for figure 1
library(tidyr)
data <- data[,c(1,2)]
data$Species <- paste0("B. ", data$Species)
datatab <- as.data.frame(t(table(data)))

orderedeffects<-as.factor(datatab$Location)[c(5,8,9,6,1,7,3,2,4)]

datatab$Location <-as.factor(datatab$Location)
datatab$Location <-factor(datatab$Location,levels=orderedeffects)

coeffplot <- ggplot(datatab, aes(fill=Species, y=Freq, x=Location)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colours) +
  scale_y_continuous(expand = c(0,0))+
  xlab("") +
  ylab("Number of individuals") +
  theme_bw()+
  theme(text = element_text(size = 15)) +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(face = "italic"),
        legend.position = "none")
##

tiff("~/Desktop/HostByLocNoLegend.tif", width = 297, height = 210, res = 300, units = "mm")
coeffplot
dev.off()

###

mat<-matrix(NA,nrow=length(unique(reshapeddata$Species))*length(unique(reshapeddata$Virus)),ncol=8)
colnames(mat)<-c("Virus","Location","Host","LL","L","M","H","HH")
mat<-as.data.frame(mat)

k<-1
for (i in 1:length(unique(reshapeddata$Virus))) {
  for (u in 1:length(unique(reshapeddata$Host))) {
    temp<-subset(reshapeddata,Virus==unique(reshapeddata$Virus)[i]&Host==unique(reshapeddata$Host)[u])
    temp<-temp[!(temp$Status%in%NA),]
    mat[k,1]<-as.character(unique(reshapeddata$Virus)[i])
    mat[k,3]<-as.character(unique(reshapeddata$Host)[u])
    if (nrow(temp)>0){
      message(nrow(temp))
      sims<-rbeta(20000,(1+sum(temp$Status)),(1+(nrow(temp)-sum(temp$Status))))
      mat[k,6]<-posterior.mode(as.mcmc(sims))
      mat[k,c(4,8)]<-SPIn(as.mcmc(sims),conf=0.9,lb=0,ub=1)$spin
      mat[k,c(5,7)]<-SPIn(as.mcmc(sims),conf=0.5,lb=0,ub=1)$spin
      print(k)
    }
    k<-k+1
  }
}  

mat
mat$Virus[mat$Virus == "SBPV"] <- "Slow bee paralysis virus"
mat$Virus[mat$Virus == "ABPV"] <- "Acute bee paralysis virus"
mat$Virus <- as.factor(mat$Virus)
mat$Virus<-factor(mat$Virus, levels = levels(mat$Virus)[c(3,5,6,4,2,1)])

library(ggplot2)

colours <- c("#68023F", "#008169", "#EF0096", "#00DCB5", "#FFCFE2", "#003C86", "#9400E6", "#009FFA",
             "#FF71FD", "#7CFFFA", "#6A0213", "#008607", "#00E307")

coeffplot <- ggplot(mat)
coeffplot <- coeffplot + facet_wrap(~Virus, ncol=3)
coeffplot <- coeffplot + geom_point(aes(x = Host, y = M, colour = Host), position = position_dodge(width = 1/2), size=2)
coeffplot <- coeffplot + geom_linerange(aes(x = Host, ymin = L,ymax = H, colour = Host),lwd = 2, position = position_dodge(width = 1/2))
coeffplot <- coeffplot + geom_linerange(aes(x = Host, ymin = LL,ymax = HH, colour = Host),lwd = 1, position = position_dodge(width = 1/2))
coeffplot <- coeffplot + theme_bw() + scale_color_manual(values = colours) + theme(legend.position = "none") + xlab("Host species")
coeffplot <- coeffplot + ylab("Prevalence") + theme(axis.text.x = element_text(angle = 45, hjust = 1,face = "italic", size = 9), 
                                                    axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
print(coeffplot)  # The trick to these is position_dodge()

tiff("~/Desktop/Host.tif", width = 267.3, height = 189, res = 300, units = "mm")
coeffplot
dev.off()