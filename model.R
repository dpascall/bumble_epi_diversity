library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
rm(list = ls())

data <- read.csv("~/Desktop/BumblebeePaper1/P1datafinal.csv", stringsAsFactors = F)[,c(2:ncol(read.csv("~/Desktop/BumblebeePaper1/P1datafinal.csv")))]
data$LocationSpecies <- paste(data$Location, data$Species, sep = ".")

speciesmap <- cbind(unique(data$Species),as.numeric(as.factor(unique(data$Species))))
locationmap <- cbind(unique(data$Location),as.numeric(as.factor(unique(data$Location))))
specieslocationmap <- cbind(unique(paste(data$Location, data$Species, sep = ".")), as.numeric(as.factor(unique(paste(data$Location, data$Species, sep = ".")))))

LocationSpeciesNumber <- NULL
SpeciesNumber <- NULL
LocationNumber <- NULL
for (i in 1:nrow(data)) {
  LocationSpeciesNumber <- as.numeric(c(LocationSpeciesNumber, specieslocationmap[which(specieslocationmap[,1] == data$LocationSpecies[i]) ,2]))
  SpeciesNumber <- as.numeric(c(SpeciesNumber, speciesmap[which(speciesmap[,1] == data$Species[i]) ,2]))
  LocationNumber <- as.numeric(c(LocationNumber, locationmap[which(locationmap[,1] == data$Location[i]) ,2]))
}
data$LocationSpeciesNumber <- LocationSpeciesNumber
data$SpeciesNumber <- SpeciesNumber
data$LocationNumber <- LocationNumber


##read in files for prediction
AllPascuorum <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/AllPascuorum.csv")[,c(2,3)]
AllTerrestris <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/AllTerrestris.csv")[,c(2,3)]
AllLucorum <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/AllLucorum.csv")[,c(2,3)]
AllOther <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/AllOther.csv")[,c(2,3)]

ABPVPascuorum <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/ABPVPascuorum.csv")[,c(2,3)]
ABPVTerrestris <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/ABPVTerrestris.csv")[,c(2,3)]
ABPVLucorum <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/ABPVLucorum.csv")[,c(2,3)]
ABPVOther <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/ABPVOther.csv")[,c(2,3)]

SBPVPascuorum <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/SBPVPascuorum.csv")[,c(2,3)]
SBPVTerrestris <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/SBPVTerrestris.csv")[,c(2,3)]
SBPVLucorum <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/SBPVLucorum.csv")[,c(2,3)]
SBPVOther <- read.csv("~/Desktop/BumblebeePaper1/PredictionTargets/SBPVOther.csv")[,c(2,3)]

##put them in list for purposes of iteration
listedpredictions <- list(AllPascuorum, AllTerrestris, AllLucorum, AllOther, ABPVPascuorum, ABPVTerrestris, ABPVLucorum, ABPVOther, SBPVPascuorum, SBPVTerrestris, SBPVLucorum, SBPVOther)

for (i in 1:length(listedpredictions)) {
  listedpredictions[[i]] <- cbind(listedpredictions[[i]], rep(NA, nrow(listedpredictions[[i]])), rep(NA, nrow(listedpredictions[[i]])), rep(NA, nrow(listedpredictions[[i]])), rep(NA, nrow(listedpredictions[[i]])), rep(NA, nrow(listedpredictions[[i]])), rep(NA, nrow(listedpredictions[[i]])))
  for (j in 1:nrow(listedpredictions[[i]])) {
    if (as.character(listedpredictions[[i]][j,1]) %in% speciesmap[,1]) {
      listedpredictions[[i]][j,3] <- speciesmap[which(speciesmap[,1] == as.character(listedpredictions[[i]][j,1])),2]
    }
    if (as.character(listedpredictions[[i]][j,2]) %in% locationmap[,1]) {
      listedpredictions[[i]][j,4] <- locationmap[which(locationmap[,1] == as.character(listedpredictions[[i]][j,2])),2]
    }
    if (paste(as.character(listedpredictions[[i]][j,2]), as.character(listedpredictions[[i]][j,1]), sep = ".") %in% specieslocationmap[,1]) {
      listedpredictions[[i]][j,5] <- specieslocationmap[which(specieslocationmap[,1] == paste(as.character(listedpredictions[[i]][j,2]), as.character(listedpredictions[[i]][j,1]), sep = ".")),2]
    }
    #get maxtemp, precipitation, wind speed
    listedpredictions[[i]][j,c(6,7,8)] <- subset(data, data$Location == as.character(listedpredictions[[i]][j,2]))[1, c(9,11,13)]
  }
}

listedpredictions[[4]][is.na(listedpredictions[[4]][,5]),5] <- max(LocationSpeciesNumber) + 1

data_reduced <- data[is.na(data$SBPV),]
data_SBPV <- data[!is.na(data$SBPV) & is.na(data$ABPV),]
data_all <- data[!is.na(data$ABPV),]

##generate model matrices - maxtemp, precipitation, wind speed
predmat_reduced <- array(NA,dim = c(nrow(data_reduced), 3))
for (i in 1:nrow(data_reduced)) {
  predmat_reduced[i,c(1:3)] <- as.numeric(data_reduced[i,c(9,11,13)])
}

predmat_SBPV <- array(NA,dim = c(nrow(data_SBPV), 3))
for (i in 1:nrow(data_SBPV)) {
  predmat_SBPV[i,c(1:3)] <- as.numeric(data_SBPV[i,c(9,11,13)])
}

predmat_all <- array(NA,dim = c(nrow(data_all), 3))
for (i in 1:nrow(data_all)) {
  predmat_all[i,c(1:3)] <- as.numeric(data_all[i,c(9,11,13)])
}

##Xs for prediction
preditionmodelmatrices <- list(length = 12)

for (i in 1:12) {
  preditionmodelmatrices[[i]] <- array(NA,dim = c(nrow(listedpredictions[[i]]), 3))
  for (j in 1:nrow(listedpredictions[[i]])) {
    preditionmodelmatrices[[i]][j, c(1:3)] <- as.numeric(listedpredictions[[i]][j,c(6:8)])
  }
}

moddat <- list(nY_all = nrow(data_all), 
               nY_SBPV = nrow(data_SBPV), 
               nY_reduced = nrow(data_reduced), 
               nV_all = 6, 
               nV_SBPV = 5, 
               nV_reduced = 4, 
               nL = length(unique(data$Location)), 
               nS = length(unique(data$Species)), 
               L_all = data_all$LocationNumber, 
               L_SBPV = data_SBPV$LocationNumber, 
               L_reduced = data_reduced$LocationNumber, 
               S_all = data_all$SpeciesNumber, 
               S_SBPV = data_SBPV$SpeciesNumber, 
               S_reduced = data_reduced$SpeciesNumber, 
               Y_all = as.matrix(data_all[,c(16:21)]),  
               Y_SBPV = as.matrix(data_SBPV[,c(16:20)]), 
               Y_reduced = as.matrix(data_reduced[,c(16:19)]), 
               nLS = length(unique(data$LocationSpecies)), 
               LS_all = data_all$LocationSpeciesNumber, 
               LS_SBPV = data_SBPV$LocationSpeciesNumber, 
               LS_reduced = data_reduced$LocationSpeciesNumber, 
               X_all = predmat_all, 
               X_SBPV = predmat_SBPV, 
               X_reduced = predmat_reduced, 
               nX = dim(predmat_all)[2],
               nY_pasc_all = nrow(listedpredictions[[1]]), 
               nY_pasc_SBPV = nrow(listedpredictions[[9]]), 
               nY_pasc_ABPV = nrow(listedpredictions[[5]]), 
               nY_luc_all = nrow(listedpredictions[[3]]), 
               nY_luc_SBPV = nrow(listedpredictions[[11]]), 
               nY_luc_ABPV = nrow(listedpredictions[[7]]),
               nY_terr_all = nrow(listedpredictions[[2]]),
               nY_terr_SBPV = nrow(listedpredictions[[10]]),
               nY_terr_ABPV = nrow(listedpredictions[[6]]),
               nY_other_all = nrow(listedpredictions[[4]]),
               nY_other_SBPV = nrow(listedpredictions[[12]]),
               nY_other_ABPV = nrow(listedpredictions[[8]]),
               X_pred_pasc_all = preditionmodelmatrices[[1]],
               X_pred_pasc_SBPV = preditionmodelmatrices[[9]],
               X_pred_pasc_ABPV = preditionmodelmatrices[[5]],
               X_pred_luc_all = preditionmodelmatrices[[3]],
               X_pred_luc_SBPV = preditionmodelmatrices[[11]],
               X_pred_luc_ABPV = preditionmodelmatrices[[7]],
               X_pred_terr_all = preditionmodelmatrices[[2]],
               X_pred_terr_SBPV = preditionmodelmatrices[[10]],
               X_pred_terr_ABPV = preditionmodelmatrices[[6]],
               X_pred_other_all = preditionmodelmatrices[[4]],
               X_pred_other_SBPV = preditionmodelmatrices[[12]],
               X_pred_other_ABPV = preditionmodelmatrices[[8]],
               L_pred_pasc_all = as.numeric(listedpredictions[[1]][,4]),
               L_pred_pasc_SBPV = as.numeric(listedpredictions[[9]][,4]),
               L_pred_pasc_APBV = as.numeric(listedpredictions[[5]][,4]),
               L_pred_luc_all = as.numeric(listedpredictions[[3]][,4]),
               L_pred_luc_SBPV = as.numeric(listedpredictions[[11]][,4]),
               L_pred_luc_APBV = as.numeric(listedpredictions[[7]][,4]),
               L_pred_terr_all = as.numeric(listedpredictions[[2]][,4]),
               L_pred_terr_SBPV = as.numeric(listedpredictions[[10]][,4]),
               L_pred_terr_APBV = as.numeric(listedpredictions[[6]][,4]),
               L_pred_other_all = as.numeric(listedpredictions[[4]][,4]),
               L_pred_other_SBPV = as.numeric(listedpredictions[[12]][,4]),
               L_pred_other_APBV = as.numeric(listedpredictions[[8]][,4]),
               S_pred_pasc_all = as.numeric(listedpredictions[[1]][,3]),
               S_pred_pasc_SBPV = as.numeric(listedpredictions[[9]][,3]),
               S_pred_pasc_APBV = as.numeric(listedpredictions[[5]][,3]),
               S_pred_luc_all = as.numeric(listedpredictions[[3]][,3]),
               S_pred_luc_SBPV = as.numeric(listedpredictions[[11]][,3]),
               S_pred_luc_APBV = as.numeric(listedpredictions[[7]][,3]),
               S_pred_terr_all = as.numeric(listedpredictions[[2]][,3]),
               S_pred_terr_SBPV = as.numeric(listedpredictions[[10]][,3]),
               S_pred_terr_APBV = as.numeric(listedpredictions[[6]][,3]),
               S_pred_other_all = as.numeric(listedpredictions[[4]][,3]),
               S_pred_other_SBPV = as.numeric(listedpredictions[[12]][,3]),
               S_pred_other_APBV = as.numeric(listedpredictions[[8]][,3]),
               LS_pred_pasc_all = as.numeric(listedpredictions[[1]][,5]),
               LS_pred_pasc_SBPV = as.numeric(listedpredictions[[9]][,5]),
               LS_pred_pasc_APBV = as.numeric(listedpredictions[[5]][,5]),
               LS_pred_luc_all = as.numeric(listedpredictions[[3]][,5]),
               LS_pred_luc_SBPV = as.numeric(listedpredictions[[11]][,5]),
               LS_pred_luc_APBV = as.numeric(listedpredictions[[7]][,5]),
               LS_pred_terr_all = as.numeric(listedpredictions[[2]][,5]),
               LS_pred_terr_SBPV = as.numeric(listedpredictions[[10]][,5]),
               LS_pred_terr_APBV = as.numeric(listedpredictions[[6]][,5]),
               LS_pred_other_all = as.numeric(listedpredictions[[4]][,5]),
               LS_pred_other_SBPV = as.numeric(listedpredictions[[12]][,5]),
               LS_pred_other_APBV = as.numeric(listedpredictions[[8]][,5]))
str(moddat)

fit <- stan(file = '~/Desktop/BumblebeePaper1/modelfinalallwithbetasnoncentered.stan', data = moddat, iter = 13000, warmup = 10000, thin = 1, chains = 4, control = list(max_treedepth = 40, adapt_delta = 0.9999999999999999, metric = "dense_e"), init_r = 0.5, verbose = T, refresh = 10)
save(fit, file = "~/Desktop/BumblebeePaper1/StanModelOutput.Rdata")

##generate HPD intervals

library(SPIn)

correlations <- as.matrix(fit)[,grep("Omega", colnames(as.matrix(fit)))]
correlations <- correlations[,grep("L_Omega", colnames(correlations), invert = T)]
correlations <- correlations[,-c(1,8,15,22,29,36)]
colMeans(correlations)

errors <- matrix(NA, nrow = 30, ncol = 2)

for (i in 1:ncol(correlations)) {
  errors[i,] <- SPIn(correlations[,i], conf = 0.90, lb = -1, ub = 1)$spin
}

cbind(colMeans(correlations), errors)

