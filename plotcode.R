library(rstan)
library(bayesplot)
rm(list = ls())

load("~/Desktop/StanModelOutput2020-12-12 05:58:00.Rdata")
data<-read.csv("~/Dropbox/P1datafinal.csv")[c(1:759),]
sims<-rstan::extract(fit,pars=c("alpha","beta","sigmaL","sigmaS", "sigmaLS","muL","muS"))
Omegasims<-extract(fit,pars=c("Omega"))

betasign<-matrix(NA,nrow=18,ncol=3)
betasign<-as.data.frame(betasign)

betas<-c("Precipitation","Maximum Temperature","Wind Speed")
viruses<-c("River Luinaeg virus","Loch Morlich virus","Mayfield virus 1","Mayfield virus 2", "Slow bee paralysis virus", "Acute bee paralysis virus")

k<-1
for (v in 1:6) {
  for (b in 1:3) {
    betasign[k,1]<-viruses[v]
    betasign[k,2]<-betas[b]
    betasign[k,3]<-(sum(sims$beta[,v,b]>0)/length(sims$beta[,v,b]))
    k<-k+1
  }
}

colnames(betasign)<-c("Virus","Effect","Deviation")

betasign$Virus<-as.factor(betasign$Virus)
betasign$Effect<-as.factor(betasign$Effect)

write.csv(betasign,"ProbabiltyofPositivity.csv")

coeffplot <- ggplot(betasign)
coeffplot <- coeffplot + facet_wrap(~Virus) + geom_hline(yintercept = 0.5, colour = gray(1/2), lty = 2) + geom_hline(yintercept = 1, colour = gray(1/2), lty = 2) + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + ylim(c(0,1))
coeffplot <- coeffplot + geom_point(aes(x = Effect,y =Deviation),position = position_dodge(width = 1/2),size=3)
coeffplot <- coeffplot + theme(text = element_text(size=20)) + theme_bw()
coeffplot <- coeffplot + ylab("Probabilty that sign is positive") + xlab("") + theme(axis.text.x = element_text(angle = 57, hjust = 1, size=8))
print(coeffplot)  # The trick to these is position_dodge()

omega<-array(NA,dim=c(4,4,3))
colnames(omega)<-rownames(omega)<-c("River Luinaeg virus","Loch Morlich virus","Mayfield virus 1","Mayfield virus 2")
for (v in 2:4) {
  for (v2 in 1:(v-1)){
    omega[v,v2,c(1,3)]<-SPIn(as.mcmc(Omegasims$Omega[,v,v2]),0.9)$spin
    omega[v,v2,2]<-posterior.mode(as.mcmc(Omegasims$Omega[,v,v2]))
  }
}

omega<-round(omega,digits=3)

write.csv(omega,"correlationsfrommodel.csv")

pars<-matrix(NA,nrow=174,ncol=8)
pars<-as.data.frame(pars)
viruses<-c("River Luinaeg virus","Loch Morlich virus","Mayfield virus 1","Mayfield virus 2", "Slow bee paralysis virus", "Acute bee paralysis virus")
k<-1
for (v in 1:6) {
  pars[k,1]<-viruses[v]
  pars[k,2]<-"Intercept"
  pars[k,c(3,7)]<-SPIn(as.mcmc(sims$alpha[,v]),0.9)$spin
  pars[k,c(4,6)]<-SPIn(as.mcmc(sims$alpha[,v]),0.5)$spin
  pars[k,5]<-posterior.mode(as.mcmc(sims$alpha[,v]))
  pars[k,8]<-"Fixed Effects"
  k<-k+1
}

betas<-c("Precipitation","Maximum Temperature","Wind Speed")

for (v in 1:6) {
  for (b in 1:3) {
    pars[k,1]<-viruses[v]
    pars[k,2]<-betas[b]
    pars[k,c(3,7)]<-SPIn(as.mcmc(sims$beta[,v,b]),0.9)$spin
    pars[k,c(4,6)]<-SPIn(as.mcmc(sims$beta[,v,b]),0.5)$spin
    pars[k,5]<-posterior.mode(as.mcmc(sims$beta[,v,b]))
    pars[k,8]<-"Fixed Effects"
    k<-k+1
  }
}

m<-as.data.frame(cbind(as.character(data$Location),as.numeric(as.factor(data$Location))))
m$V1<-as.character(m$V1)
m$V2<-as.character(m$V2)
m<-unique(m)
m<-m[order(as.numeric(m$V2)),]
loc<-as.character(m$V1)

for (v in 1:6) {
  for (l in 1:9) {
    pars[k,1]<-viruses[v]
    pars[k,2]<-loc[l]
    pars[k,c(3,7)]<-SPIn(as.mcmc(sims$muL[,l,v]),0.9)$spin
    pars[k,c(4,6)]<-SPIn(as.mcmc(sims$muL[,l,v]),0.5)$spin
    pars[k,5]<-posterior.mode(as.mcmc(sims$muL[,l,v]))
    pars[k,8]<-"Partially Pooled Location Effects"
    k<-k+1
  }
}

m<-as.data.frame(cbind(as.character(data$Species),as.numeric(as.factor(data$Species))))
m$V1<-as.character(m$V1)
m$V2<-as.character(m$V2)
m<-unique(m)
m<-m[order(as.numeric(m$V2)),]
spec<-paste("B.",as.character(m$V1),sep=" ")

for (v in 1:6) {
  for (s in 1:13) {
    pars[k,1]<-viruses[v]
    pars[k,2]<-spec[s]
    pars[k,c(3,7)]<-SPIn(as.mcmc(sims$muS[,s,v]),0.9)$spin
    pars[k,c(4,6)]<-SPIn(as.mcmc(sims$muS[,s,v]),0.5)$spin
    pars[k,5]<-posterior.mode(as.mcmc(sims$muS[,s,v]))
    pars[k,8]<-"Partially Pooled Species Effects"
    k<-k+1
  }
}

for (v in 1:6) {
  pars[k,1]<-viruses[v]
  pars[k,2]<-"Species"
  pars[k,c(3,7)]<-SPIn(as.mcmc(sims$sigmaS[,v]),0.9)$spin
  pars[k,c(4,6)]<-SPIn(as.mcmc(sims$sigmaS[,v]),0.5)$spin
  pars[k,5]<-posterior.mode(as.mcmc(sims$sigmaS[,v]))
  pars[k,8]<-"Variances"
  k<-k+1
}

for (v in 1:6) {
  pars[k,1]<-viruses[v]
  pars[k,2]<-"Location"
  pars[k,c(3,7)]<-SPIn(as.mcmc(sims$sigmaL[,v]),0.9)$spin
  pars[k,c(4,6)]<-SPIn(as.mcmc(sims$sigmaL[,v]),0.5)$spin
  pars[k,5]<-posterior.mode(as.mcmc(sims$sigmaL[,v]))
  pars[k,8]<-"Variances"
  k<-k+1
}

for (v in 1:6) {
  pars[k,1]<-viruses[v]
  pars[k,2]<-"Species-Location Interaction"
  pars[k,c(3,7)]<-SPIn(as.mcmc(sims$sigmaLS[,v]),0.9)$spin
  pars[k,c(4,6)]<-SPIn(as.mcmc(sims$sigmaLS[,v]),0.5)$spin
  pars[k,5]<-posterior.mode(as.mcmc(sims$sigmaLS[,v]))
  pars[k,8]<-"Variances"
  k<-k+1
}

colnames(pars)<-c("Virus","Parameter","LL","L","M","H","HH","Type")

orderedeffects<-pars[c(169,163,158,91:79,33:25,9,8,7,1),2]

pars$Parameter<-as.factor(pars$Parameter)
pars$Parameter<-factor(pars$Parameter,levels=orderedeffects)

pars$Virus <- as.factor(pars$Virus)
pars$Virus<-factor(pars$Virus,levels=levels(pars$Virus)[c(3,4,5,2,6,1)])

pars <- pars[pars$Type %in% c("Fixed Effects", "Variances"),]

coeffplot <- ggplot(pars)
coeffplot <- coeffplot + facet_grid(Type~Virus,scales = "free",space = "free_y") + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) #+ ylim(c(-7,7))
coeffplot <- coeffplot + geom_point(aes(x = Parameter, y = M),size = 2.5, position = position_dodge(width = 1/2))
coeffplot <- coeffplot + geom_linerange(aes(x = Parameter,ymin = LL,ymax = HH),lwd = 1, position = position_dodge(width = 1/2))
coeffplot <- coeffplot + geom_linerange(aes(x = Parameter, ymin = L,ymax = H),lwd = 2, position = position_dodge(width = 1/2))
coeffplot <- coeffplot  + theme_bw() + coord_flip() + ylim(c(-4.5,4.5))
coeffplot <- coeffplot + ylab("Effect Size") + xlab("") + theme(text = element_text(size=11))
print(coeffplot)  # The trick to these is position_dodge()

png("~/Desktop/ModelOutputNoIt.png", width = 297, height = 100, res = 600, units = "mm")
coeffplot
dev.off()
