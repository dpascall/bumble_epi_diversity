library(gtools)
library(SPIn)
library(MGLM)
library(rdiversity)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(vegan)
rm(list=ls())

data <- read.csv("~/Desktop/bumble_epi_diversity/finaldata.csv")[c(1:759),]
tabcommunity <- table(data$Location,data$Species)

##add muscorum

tabcommunity <- cbind(tabcommunity, c(0,0,0,0,29,0,0,30,0))

reps <- array(NA,dim=c(ncol(tabcommunity)+10,nrow(tabcommunity),10000))
distances <- array(NA, dim=c(9,9,10000))

k<-1
for (i in 1:9) {
    reps[,i,] <- t(rdirichlet(10000,c(as.vector(tabcommunity[i,]+1), rep(1,10))))
    print(k)
    k<-k+1
}

for (i in 1:10000) {
	distances[,,i] <- as.matrix(vegdist(t(reps[,,i]), method = "horn"))
}

meandistances <- apply(distances, c(1,2), mean)
colnames(meandistances) <- row.names(meandistances) <- row.names(tabcommunity)

summaryindexes <- array(NA,dim=c(nrow(tabcommunity),nrow(tabcommunity),3))

for (i in 1:8) {
  for (q in (i+1):9) {
    summaryindexes[q,i,c(1,3)]<-SPIn(distances[q,i,],conf=0.9,lb=0,ub=1)$spin
  }
}

summaryindexes[,,2] <- apply(distances, c(1,2), mean)

diag(summaryindexes[,,1])<-0
diag(summaryindexes[,,2])<-0
diag(summaryindexes[,,3])<-0

colnames(summaryindexes)<-rownames(summaryindexes)<-rownames(tabcommunity)
summaryindexes<-round(summaryindexes,digits=3)
write.csv(summaryindexes[,,],file="~/Dropbox/Dissimilaritywithuncertainty.csv")
