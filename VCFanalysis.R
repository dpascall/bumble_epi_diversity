library(vcfR)
rm(list=ls())

bomb <- read.vcfR("~/Desktop/David/Desktop/Pre-Glasgow/Raw/Merged/varsbomb.vcf")@fix
luc <-  read.vcfR("~/Desktop/David/Desktop/Pre-Glasgow/Raw/Merged/varsluc.vcf")@fix
pasc <-  read.vcfR("~/Desktop/David/Desktop/Pre-Glasgow/Raw/Merged/varspasc.vcf")@fix
terr <-  read.vcfR("~/Desktop/David/Desktop/Pre-Glasgow/Raw/Merged/varster.vcf")@fix

bombdepth <- do.call("c", lapply(lapply(bomb[,8], function (x) {strsplit(x,"[;]")[[1]][[1]]}), function(x){strsplit(x,"[=]")[[1]][[2]]}))
lucdepth <- do.call("c", lapply(lapply(luc[,8], function (x) {strsplit(x,"[;]")[[1]][[1]]}), function(x){strsplit(x,"[=]")[[1]][[2]]}))
terrdepth <- do.call("c", lapply(lapply(terr[,8], function (x) {strsplit(x,"[;]")[[1]][[1]]}), function(x){strsplit(x,"[=]")[[1]][[2]]}))
pascdepth <- do.call("c", lapply(lapply(pasc[,8], function (x) {strsplit(x,"[;]")[[1]][[1]]}), function(x){strsplit(x,"[=]")[[1]][[2]]}))

bomb <- as.data.frame(bomb)
luc <- as.data.frame(luc)
terr <- as.data.frame(terr)
pasc <- as.data.frame(pasc)

bomb$depth <- as.numeric(bombdepth)
luc$depth <- as.numeric(lucdepth)
terr$depth <- as.numeric(terrdepth)
pasc$depth <- as.numeric(pascdepth)

summaryb <- bomb %>% dplyr::group_by(CHROM) %>%
  dplyr::summarise(median <- median(depth),
            lower <- range(depth)[1],
            upper <- range(depth)[2]) %>%
  data.frame();summaryb

summaryt <- terr %>% dplyr::group_by(CHROM) %>%
  dplyr::summarise(median <- median(depth),
                   lower <- range(depth)[1],
                   upper <- range(depth)[2]) %>%
  data.frame();summaryt

summaryl <- luc %>% dplyr::group_by(CHROM) %>%
  dplyr::summarise(median <- median(depth),
                   lower <- range(depth)[1],
                   upper <- range(depth)[2]) %>%
  data.frame();summaryl

summaryp <- pasc %>% dplyr::group_by(CHROM) %>%
  dplyr::summarise(median <- median(depth),
                   lower <- range(depth)[1],
                   upper <- range(depth)[2]) %>%
  data.frame();summaryp

bomb$AF <- as.numeric(do.call("c", lapply(lapply(as.character(bomb[,8]), function (x) {strsplit(x,"[;]")[[1]][[2]]}), function(x){strsplit(x,"[=]")[[1]][[2]]})))
bomb <- bomb[!bomb$AF == 1,]
bomb <- bomb[as.numeric(as.character(bomb$POS))>=200,]
bombLMV <- bomb[as.numeric(as.character(bomb$POS))<=1301 & bomb$CHROM == "Loch_Morlich_virus",]
bombRLV <- bomb[as.numeric(as.character(bomb$POS))<=1336 & bomb$CHROM == "River_Liunaeg_virus",]
bombMV2 <- bomb[as.numeric(as.character(bomb$POS))<=1283 & bomb$CHROM == "Mayfield_virus_2",]
bombMV1 <- bomb[as.numeric(as.character(bomb$POS))<=1283 & bomb$CHROM == "Mayfield_virus_1",]
bombSBPV <- bomb[as.numeric(as.character(bomb$POS))<=1319 & bomb$CHROM == "EU035616.1_Slow_bee_paralysis_virus_strain_Rothamsted_polyprotein_gene,_complete_cds",]
bombABPV <- bomb[as.numeric(as.character(bomb$POS))<=1322 & bomb$CHROM == "AF150629.1_Acute_bee_paralysis_virus,_complete_genome,_complete_cds",]

94/(1301-200)*100 #LMV
0 #ABPV
43/(1283-200)*100 #MV1
126/(1283-200)*100 #MV2
59/(1319-200)*100 #SBPV
125/(1336-200)*100 #RLV

luc$AF <- as.numeric(do.call("c", lapply(lapply(as.character(luc[,8]), function (x) {strsplit(x,"[;]")[[1]][[2]]}), function(x){strsplit(x,"[=]")[[1]][[2]]})))
luc <- luc[!luc$AF == 1,]
luc <- luc[as.numeric(as.character(luc$POS))>=200,]
lucMV1 <- luc[as.numeric(as.character(luc$POS))<=1283 & luc$CHROM == "Mayfield_virus_1",]
lucSBPV <- luc[as.numeric(as.character(luc$POS))<=1319 & luc$CHROM == "EU035616.1_Slow_bee_paralysis_virus_strain_Rothamsted_polyprotein_gene,_complete_cds",]
lucABPV <- luc[as.numeric(as.character(luc$POS))<=1322 & luc$CHROM == "AF150629.1_Acute_bee_paralysis_virus,_complete_genome,_complete_cds",]

0 #ABPV
93/(1283-200)*100 #MV1
17/(1319-200)*100 #SBPV

terr$AF <- as.numeric(do.call("c", lapply(lapply(as.character(terr[,8]), function (x) {strsplit(x,"[;]")[[1]][[2]]}), function(x){strsplit(x,"[=]")[[1]][[2]]})))
terr <- terr[!terr$AF == 1,]
terr <- terr[as.numeric(as.character(terr$POS))>=200,]
terrMV2 <- terr[as.numeric(as.character(terr$POS))<=1283 & terr$CHROM == "Mayfield_virus_2",]
terrMV1 <- terr[as.numeric(as.character(terr$POS))<=1283 & terr$CHROM == "Mayfield_virus_1",]
terrSBPV <- terr[as.numeric(as.character(terr$POS))<=1319 & terr$CHROM == "EU035616.1_Slow_bee_paralysis_virus_strain_Rothamsted_polyprotein_gene,_complete_cds",]

73/(1283-200)*100 #MV1
12/(1319-200)*100 #SBPV
44/(1283-200)*100 #MV2


pasc$AF <- as.numeric(do.call("c", lapply(lapply(as.character(pasc[,8]), function (x) {strsplit(x,"[;]")[[1]][[2]]}), function(x){strsplit(x,"[=]")[[1]][[2]]})))
pasc <- pasc[!pasc$AF == 1,]
pasc <- pasc[as.numeric(as.character(pasc$POS))>=200,]
pascSBPV <- pasc[as.numeric(as.character(pasc$POS))<=1319 & pasc$CHROM == "EU035616.1_Slow_bee_paralysis_virus_strain_Rothamsted_polyprotein_gene,_complete_cds",]
pascABPV <- pasc[as.numeric(as.character(pasc$POS))<=1322 & pasc$CHROM == "AF150629.1_Acute_bee_paralysis_virus,_complete_genome,_complete_cds",]
pascMV2 <- pasc[as.numeric(as.character(pasc$POS))<=1283 & pasc$CHROM == "Mayfield_virus_2",]

0 #ABPV
42/(1319-200)*100 #SBPV
128/(1283-200)*100 #MV2

rm("bomb", "luc", "bombdepth", "lucdepth", "pasc", "pascdepth", "summaryb", "summaryl", "summaryp", "summaryt", "terr", "terrdepth")

watthetaraw <- matrix(NA,ncol = 4, nrow = 6)

data <- read.csv("~/Dropbox/P1datafinal.csv", stringsAsFactors = F)[,c(2:ncol(read.csv("~/Dropbox/P1datafinal.csv")))]

tdata <- data[!data$BGI.pot == "" & data$Species == "terrestris", c("RLV","LMV","MV1","MV2","SBPV","ABPV")]
pdata <- data[!data$BGI.pot == "" & data$Species == "pascuorum", c("RLV","LMV","MV1","MV2","SBPV","ABPV")]
ldata <- data[!data$BGI.pot == "" & data$Species == "lucorum", c("RLV","LMV","MV1","MV2","SBPV","ABPV")]
odata <- data[!data$BGI.pot == "" & !data$Species %in% c("lucorum", "pascuorum", "terrestris"), c("RLV","LMV","MV1","MV2","SBPV","ABPV")]

library(rstan)
load("~/Desktop/StanModelOutput2020-12-12 05:58:00.Rdata")

extras <- as.matrix(fit)[,colnames(as.matrix(fit))%in%c("pasc_RLV", "pasc_LMV", "pasc_MV1", "pasc_MV2", "pasc_SBPV" , "pasc_ABPV",
                                "terr_RLV", "terr_LMV", "terr_MV1", "terr_MV2", "terr_SBPV", "terr_ABPV",
                                "luc_RLV", "luc_LMV", "luc_MV1", "luc_MV2", "luc_SBPV", "luc_ABPV",
                                "other_RLV", "other_LMV", "other_MV1", "other_MV2", "other_SBPV", "other_ABPV")]

bounds <- matrix(NA, nrow = ncol(extras), ncol = 3)
row.names(bounds) <- colnames(extras)
bounds[,1] <- apply(extras, 2, median)
for (i in 1:nrow(bounds)) {
  bounds[i, c(2,3)] <- quantile(extras[,i], probs = c(0.05, 0.95))
}

tsums_med <- as.numeric(colSums(tdata, na.rm = T)) + bounds[c(7:12), 1]
psums_med <- as.numeric(colSums(pdata, na.rm = T)) + bounds[c(1:6), 1]
lsums_med <- as.numeric(colSums(ldata, na.rm = T)) + bounds[c(13:18), 1]
osums_med <- as.numeric(colSums(odata, na.rm = T)) + bounds[c(19:24), 1]

tsums_low <- as.numeric(colSums(tdata, na.rm = T)) + bounds[c(7:12), 2]
psums_low <- as.numeric(colSums(pdata, na.rm = T)) + bounds[c(1:6), 2]
lsums_low <- as.numeric(colSums(ldata, na.rm = T)) + bounds[c(13:18), 2]
osums_low <- as.numeric(colSums(odata, na.rm = T)) + bounds[c(19:24), 2]

tsums_upp <- as.numeric(colSums(tdata, na.rm = T)) + bounds[c(7:12), 3]
psums_upp <- as.numeric(colSums(pdata, na.rm = T)) + bounds[c(1:6), 3]
lsums_upp <- as.numeric(colSums(ldata, na.rm = T)) + bounds[c(13:18), 3]
osums_upp <- as.numeric(colSums(odata, na.rm = T)) + bounds[c(19:24), 3]

f <- function(n) {
  sum(1 / (1:(n-1)))
}

##other
nrow(bombRLV)/(f(osums_med[1])*(1501-400)) ##RLV
nrow(bombLMV)/(f(osums_med[2])*(1501-400)) ##LMV
nrow(bombMV1)/(f(osums_med[3])*(1483-400)) ##MV1
nrow(bombMV2)/(f(osums_med[4])*(1483-400)) ##MV2
nrow(bombSBPV)/(f(osums_med[5]+110)*(1519-400)) ##SBPV
nrow(bombABPV)/(f(osums_med[6])*(1522-400)) ##ABPV

nrow(bombRLV)/(f(osums_low[1])*(1501-400)) ##RLV
nrow(bombLMV)/(f(osums_low[2])*(1501-400)) ##LMV
nrow(bombMV1)/(f(osums_low[3])*(1483-400)) ##MV1
nrow(bombMV2)/(f(osums_low[4])*(1483-400)) ##MV2
nrow(bombSBPV)/(f(osums_low[5]+110)*(1519-400)) ##SBPV
nrow(bombABPV)/(f(osums_low[6])*(1522-400)) ##ABPV

nrow(bombRLV)/(f(osums_upp[1])*(1501-400)) ##RLV
nrow(bombLMV)/(f(osums_upp[2])*(1501-400)) ##LMV
nrow(bombMV1)/(f(osums_upp[3])*(1483-400)) ##MV1
nrow(bombMV2)/(f(osums_upp[4])*(1483-400)) ##MV2
nrow(bombSBPV)/(f(osums_upp[5]+110)*(1519-400)) ##SBPV
nrow(bombABPV)/(f(osums_upp[6])*(1522-400)) ##ABPV

##terr
nrow(terrMV1)/(f(tsums_med[3])*(1483-400)) ##MV1
nrow(terrMV2)/(f(tsums_med[4])*(1483-400)) ##MV2
nrow(terrSBPV)/(f(tsums_med[5])*(1519-400)) ##SBPV

nrow(terrMV1)/(f(tsums_low[3])*(1483-400)) ##MV1
nrow(terrMV2)/(f(tsums_low[4])*(1483-400)) ##MV2
nrow(terrSBPV)/(f(tsums_low[5])*(1519-400)) ##SBPV

nrow(terrMV1)/(f(tsums_upp[3])*(1483-400)) ##MV1
nrow(terrMV2)/(f(tsums_upp[4])*(1483-400)) ##MV2
nrow(terrSBPV)/(f(tsums_upp[5])*(1519-400)) ##SBPV

##luc
nrow(lucMV1)/(f(lsums_med[3])*(1483-400)) ##MV1
nrow(lucSBPV)/(f(lsums_med[5])*(1519-400)) ##SBPV
nrow(lucABPV)/(f(lsums_med[6])*(1522-400)) ##ABPV

nrow(lucMV1)/(f(lsums_low[3])*(1483-400)) ##MV1
nrow(lucSBPV)/(f(lsums_low[5])*(1519-400)) ##SBPV
nrow(lucABPV)/(f(lsums_low[6])*(1522-400)) ##ABPV

nrow(lucMV1)/(f(lsums_upp[3])*(1483-400)) ##MV1
nrow(lucSBPV)/(f(lsums_upp[5])*(1519-400)) ##SBPV
nrow(lucABPV)/(f(lsums_upp[6])*(1522-400)) ##ABPV

##pasc
nrow(pascMV2)/(f(psums_med[4])*(1483-400)) ##MV2
nrow(pascSBPV)/(f(psums_med[5])*(1519-400)) ##SBPV
nrow(pascABPV)/(f(psums_med[6])*(1522-400)) ##ABPV

nrow(pascMV2)/(f(psums_low[4])*(1483-400)) ##MV2
nrow(pascSBPV)/(f(psums_low[5])*(1519-400)) ##SBPV
nrow(pascABPV)/(f(psums_low[6])*(1522-400)) ##ABPV

nrow(pascMV2)/(f(psums_upp[4])*(1483-400)) ##MV2
nrow(pascSBPV)/(f(psums_upp[5])*(1519-400)) ##SBPV
nrow(pascABPV)/(f(psums_upp[6])*(1522-400)) ##ABPV

