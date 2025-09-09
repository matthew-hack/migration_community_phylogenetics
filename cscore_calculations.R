#cscore_calculations.R



library(ebirdst)
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(rnaturalearth)
library(geodata)
library(ggplot2)
library(fields)
library(stringr)
library(TreeTools)
library(viridis)
#library(rgdal)
library(raster)
library(devtools)
#devtools::install_github("GotelliLab/EcoSimR")
library(EcoSimR)
library(ecospat)

bcr_12_week_26 <- read.csv("output/bcr_community_data/26/community_26_bcr_12.csv")
bcr_12_week_26_presence <- bcr_12_week_26 %>% mutate_if(is.numeric, ~1 * (. != 0))
bcr_12_week_26_presence_subsample <- bcr_12_week_26_presence[sample(nrow(bcr_12_week_26_presence), nrow(bcr_12_week_26_presence)/10, replace=FALSE),]

#cscore_results <- cooc_null_model(bcr_12_week_26_presence_subsample, algo = "sim9", metric = "c_score", nReps = 100)
#plot(cscore_results,type="cooc")
#write.csv(cscore_results, "output/tables/c_scores_bcr_12_week_26.csv")


ecospot <- ecospat.Cscore(bcr_12_week_26_presence_subsample, 100, outpath="output/tables")
c_scores <- read.delim("output/tables/Cscores.txt")

c_scores[sapply(c_scores, is.infinite)] <- NA
mean(c_scores$SES_Cscore, na.rm=TRUE) #mean c score: 0.2231 (overdispersion)


#### deal with phylogenies
passerine.picCleanTree.bcr12week26 <- match.phylo.comm(phy = passerine.tree, 
                                                       comm = bcr_12_week_26_presence_subsample)$phy
passerine.picCleanComm.bcr12week26 <- match.phylo.comm(phy = passerine.tree, 
                                                       comm = bcr_12_week_26_presence_subsample)$comm

# Create a Phylogenetic Distance Matrix
passerine.cophenDist <- cophenetic.phylo(passerine.picCleanTree.bcr12week26)

c_scores$phylogenetic_distance <- rep(NA, nrow(c_scores))
for(i in 1:nrow(c_scores)){
  c_scores[i, "phylogenetic_distance"] <- passerine.cophenDist[c_scores[i,"Sp1"],c_scores[i,"Sp2"]]
}

ggplot(c_scores, aes(phylogenetic_distance, SES_Cscore)) +
  geom_point() + 
  geom_smooth(method = "loess")

mean(c_scores$phylogenetic_distance) #mean phylogenetic distance: 80.98
nrow(c_scores[c_scores$phylogenetic_distance < 50,])/nrow(c_scores) #0.272 (more closer relatives)


warblers <- c("Seiurus_aurocapilla","Seiurus_motacilla","Seiurus_noveboracensis","Vermivora_chrysoptera",
              "Mniotilta_varia","Vermivora_peregrina","Vermivora_celata","Vermivora_ruficapilla",
              "Oporornis_agilis","Oporornis_philadelphia","Geothlypis_trichas","Wilsonia_citrina",
              "Setophaga_ruticilla","Dendroica_kirtlandii","Dendroica_tigrina","Dendroica_cerulea",
              "Parula_americana","Dendroica_magnolia","Dendroica_castanea","Dendroica_fusca",
              "Dendroica_petechia","Dendroica_pensylvanica","Dendroica_striata","Dendroica_caerulescens",
              "Dendroica_palmarum","Dendroica_pinus","Dendroica_coronata","Dendroica_discolor",
              "Dendroica_virens","Wilsonia_canadensis","Wilsonia_pusilla")

warbler_c_scores <- filter(c_scores, Sp1 %in% warblers, Sp2 %in% warblers)
View(warbler_c_scores)

mean(warbler_c_scores$SES_Cscore, na.rm=TRUE) #mean warbler c score: -4.6795 (clustering)
#median warbler c score: -2.2052 

setophaga <- c("Wilsonia_citrina","Setophaga_ruticilla","Dendroica_kirtlandii","Dendroica_tigrina","Dendroica_cerulea",
               "Parula_americana","Dendroica_magnolia","Dendroica_castanea","Dendroica_fusca",
               "Dendroica_petechia","Dendroica_pensylvanica","Dendroica_striata","Dendroica_caerulescens",
               "Dendroica_palmarum","Dendroica_pinus","Dendroica_coronata","Dendroica_discolor",
               "Dendroica_virens")
setophaga_c_scores <- filter(c_scores, Sp1 %in% setophaga, Sp2 %in% setophaga)

mean(setophaga_c_scores$SES_Cscore, na.rm=TRUE) #mean setophaga c score: -3.8391 (slightly less clustered)
#median setophaga c score: -1.6407

ggplot(warbler_c_scores, aes(phylogenetic_distance, SES_Cscore)) +
  geom_point() + 
  geom_smooth(method = "loess")

ggplot(setophaga_c_scores, aes(phylogenetic_distance, SES_Cscore)) +
  geom_point() + 
  geom_smooth(method = "loess")

summary(lm(SES_Cscore ~ phylogenetic_distance, data = setophaga_c_scores))


### repeat with week 1, bcr 12

bcr_12_week_1 <- read.csv("output/bcr_community_data/1/community_1_bcr_12.csv")
bcr_12_week_1_presence <- bcr_12_week_1 %>% mutate_if(is.numeric, ~1 * (. != 0))
bcr_12_week_1_presence_subsample <- bcr_12_week_1_presence[sample(nrow(bcr_12_week_1_presence), nrow(bcr_12_week_1_presence)/100, replace=FALSE),]


ecospat_bcr_12_week_1 <- ecospat.Cscore(bcr_12_week_1_presence_subsample, 100, outpath="output/tables")
c_scores_bcr_12_week_1 <- read.delim("output/tables/Cscores.txt")

c_scores_bcr_12_week_1[sapply(c_scores_bcr_12_week_1, is.infinite)] <- NA
mean(c_scores_bcr_12_week_1$SES_Cscore, na.rm=TRUE) #mean c score: 0.382115 (overdispersion, slightly less)

passerine.picCleanTree.bcr12week1 <- match.phylo.comm(phy = passerine.tree, 
                                                      comm = bcr_12_week_1_presence_subsample)$phy
passerine.picCleanComm.bcr12week1 <- match.phylo.comm(phy = passerine.tree, 
                                                      comm = bcr_12_week_1_presence_subsample)$comm

# Create a Phylogenetic Distance Matrix
passerine.cophenDist <- cophenetic.phylo(passerine.picCleanTree.bcr12week1)

c_scores_bcr_12_week_1$phylogenetic_distance <- rep(NA, nrow(c_scores_bcr_12_week_1))
for(i in 1:nrow(c_scores_bcr_12_week_1)){
  c_scores_bcr_12_week_1[i, "phylogenetic_distance"] <- passerine.cophenDist[c_scores_bcr_12_week_1[i,"Sp1"],c_scores_bcr_12_week_1[i,"Sp2"]]
}

ggplot(c_scores_bcr_12_week_1, aes(phylogenetic_distance, SES_Cscore)) +
  geom_point() + 
  geom_smooth(method = "loess")

mean(c_scores_bcr_12_week_1$phylogenetic_distance) #mean phylogenetic distance: 77.06
nrow(c_scores_bcr_12_week_1[c_scores_bcr_12_week_1$phylogenetic_distance < 50,])/nrow(c_scores_bcr_12_week_1) #0.183 (fewer closer relatives)


### repeat in BCR 27


bcr_27_week_26 <- read.csv("output/bcr_community_data/26/community_26_bcr_27.csv")
bcr_27_week_26_presence <- bcr_27_week_26 %>% mutate_if(is.numeric, ~1 * (. != 0))
bcr_27_week_26_presence_subsample <- bcr_27_week_26_presence[sample(nrow(bcr_27_week_26_presence), nrow(bcr_27_week_26_presence)/100, replace=FALSE),]

#cscore_results <- cooc_null_model(bcr_27_week_26_presence_subsample, algo = "sim9", metric = "c_score", nReps = 100)
#plot(cscore_results,type="cooc")
#write.csv(cscore_results, "output/tables/c_scores_bcr_27_week_26.csv")


ecospot <- ecospat.Cscore(bcr_27_week_26_presence_subsample, 100, outpath="output/tables")
c_scores_bcr_27_week_26 <- read.delim("output/tables/Cscores.txt")

c_scores_bcr_27_week_26[sapply(c_scores_bcr_27_week_26, is.infinite)] <- NA
mean(c_scores_bcr_27_week_26$SES_Cscore, na.rm=TRUE) #mean c score: 0.2334 (overdispersion)

passerine.picCleanTree.bcr27week26 <- match.phylo.comm(phy = passerine.tree, 
                                                       comm = bcr_27_week_26_presence_subsample)$phy
passerine.picCleanComm.bcr27week26 <- match.phylo.comm(phy = passerine.tree, 
                                                       comm = bcr_27_week_26_presence_subsample)$comm

# Create a Phylogenetic Distance Matrix
passerine.cophenDist <- cophenetic.phylo(passerine.picCleanTree.bcr27week26)

c_scores_bcr_27_week_26$phylogenetic_distance <- rep(NA, nrow(c_scores_bcr_27_week_26))
for(i in 1:nrow(c_scores_bcr_27_week_26)){
  c_scores_bcr_27_week_26[i, "phylogenetic_distance"] <- passerine.cophenDist[c_scores_bcr_27_week_26[i,"Sp1"],c_scores_bcr_27_week_26[i,"Sp2"]]
}

ggplot(c_scores_bcr_27_week_26, aes(phylogenetic_distance, SES_Cscore)) +
  geom_point() + 
  geom_smooth(method = "loess")

mean(c_scores_bcr_27_week_26$phylogenetic_distance) #mean phylogenetic distance: 86.87
nrow(c_scores_bcr_27_week_26[c_scores_bcr_27_week_26$phylogenetic_distance < 50,])/nrow(c_scores_bcr_27_week_26) #0.243 (somewhat more closer relatives)


### repeat with week 1, bcr 27
bcr_27_week_1 <- read.csv("output/bcr_community_data/1/community_1_bcr_27.csv")
bcr_27_week_1_presence <- bcr_27_week_1 %>% mutate_if(is.numeric, ~1 * (. != 0))
bcr_27_week_1_presence_subsample <- bcr_27_week_1_presence[sample(nrow(bcr_27_week_1_presence), nrow(bcr_27_week_1_presence)/100, replace=FALSE),]

#cscore_results <- cooc_null_model(bcr_27_week_1_presence_subsample, algo = "sim9", metric = "c_score", nReps = 100)
#plot(cscore_results,type="cooc")
#write.csv(cscore_results, "output/tables/c_scores_bcr_27_week_1.csv")


ecospot <- ecospat.Cscore(bcr_27_week_1_presence_subsample, 100, outpath="output/tables")
c_scores_bcr_27_week_1 <- read.delim("output/tables/Cscores.txt")

c_scores_bcr_27_week_1[sapply(c_scores_bcr_27_week_1, is.infinite)] <- NA
mean(c_scores_bcr_27_week_1$SES_Cscore, na.rm=TRUE) #mean c score: 0.1614 (overdispersion)

passerine.picCleanTree.bcr27week1 <- match.phylo.comm(phy = passerine.tree, 
                                                      comm = bcr_27_week_1_presence_subsample)$phy
passerine.picCleanComm.bcr27week1 <- match.phylo.comm(phy = passerine.tree, 
                                                      comm = bcr_27_week_1_presence_subsample)$comm

# Create a Phylogenetic Distance Matrix
passerine.cophenDist <- cophenetic.phylo(passerine.picCleanTree.bcr27week1)

c_scores_bcr_27_week_1$phylogenetic_distance <- rep(NA, nrow(c_scores_bcr_27_week_1))
for(i in 1:nrow(c_scores_bcr_27_week_1)){
  c_scores_bcr_27_week_1[i, "phylogenetic_distance"] <- passerine.cophenDist[c_scores_bcr_27_week_1[i,"Sp1"],c_scores_bcr_27_week_1[i,"Sp2"]]
}

ggplot(c_scores_bcr_27_week_1, aes(phylogenetic_distance, SES_Cscore)) +
  geom_point() + 
  geom_smooth(method = "loess")

mean(c_scores_bcr_27_week_1$phylogenetic_distance) #mean phylogenetic distance: 77.379
nrow(c_scores_bcr_27_week_1[c_scores_bcr_27_week_1$phylogenetic_distance < 50,])/nrow(c_scores_bcr_27_week_1) #0.308 (closer relatives)





#### try to get values for mean dispersion for NRI & NTI

#BCR 12 week 26
hist(c_scores$phylogenetic_distance)
mean(c_scores$phylogenetic_distance) #80.983
sd(c_scores$phylogenetic_distance) #37.03

#calculate nearest phylogenetic neighbor for each
min(c_scores[c_scores$Sp1 == "Setophaga_ruticilla" | c_scores$Sp2 == "Setophaga_ruticilla",
             "phylogenetic_distance"])
bcr12week26_species <- colnames(bcr_12_week_26)
bcr12week26_min_phy_distances <- rep(NA, length(bcr12week26_species)) 

for(i in 1:length(bcr12week26_species)){
  bcr12week26_min_phy_distances[i] <- min(c_scores[c_scores$Sp1 == bcr12week26_species[i] |
                                                     c_scores$Sp2 == bcr12week26_species[i], 
                                                   "phylogenetic_distance"])
}

hist(bcr12week26_min_phy_distances)
mean(bcr12week26_min_phy_distances) #18.196
sd(bcr12week26_min_phy_distances) #15.432

hist(bcr12week26_min_phy_distances, freq=FALSE, col=rgb(0,0,1,1/2), 
     xlim=c(0, max(c_scores$phylogenetic_distance)), main="Phylogenetic Distances for NTI & NRI, BCR 12, Week 26",
     xlab="Phylogenetic Distance", breaks=20)
hist(c_scores$phylogenetic_distance, freq=FALSE, col=rgb(1,0,0,1/2),breaks=20,add=TRUE)
abline(v=mean(c_scores$phylogenetic_distance), col="red", lwd=3)
abline(v=mean(bcr12week26_min_phy_distances), col="blue", lwd=3)


#BCR 12 week 1
hist(c_scores_bcr_12_week_1$phylogenetic_distance)
mean(c_scores_bcr_12_week_1$phylogenetic_distance) #77.062
sd(c_scores_bcr_12_week_1$phylogenetic_distance) #27.243

#calculate nearest phylogenetic neighbor for each
bcr12week1_species <- colnames(bcr_12_week_1)
bcr12week1_min_phy_distances <- rep(NA, length(bcr12week1_species)) 

for(i in 1:length(bcr12week1_species)){
  bcr12week1_min_phy_distances[i] <- min(c_scores_bcr_12_week_1[c_scores_bcr_12_week_1$Sp1 == bcr12week1_species[i] |
                                                                  c_scores_bcr_12_week_1$Sp2 == bcr12week1_species[i], 
                                                                "phylogenetic_distance"])
}

hist(bcr12week1_min_phy_distances)
mean(bcr12week1_min_phy_distances) #25.109
sd(bcr12week1_min_phy_distances) #23.700

hist(c_scores_bcr_12_week_1$phylogenetic_distance, freq=FALSE, 
     breaks=20, col=rgb(1,0,0,1/2), main="Phylogenetic Distances for NTI & NRI, BCR 12, Week 1",
     xlab="Phylogenetic Distance")
hist(bcr12week1_min_phy_distances, freq=FALSE, col=rgb(0,0,1,1/2), 
     xlim=c(0, max(c_scores_bcr_12_week_1$phylogenetic_distance)), breaks=20, add=TRUE)
abline(v=mean(c_scores_bcr_12_week_1$phylogenetic_distance), col="red", lwd=3)
abline(v=mean(bcr12week1_min_phy_distances), col="blue", lwd=3)


#BCR 27 week 26
hist(c_scores_bcr_27_week_26$phylogenetic_distance)
mean(c_scores_bcr_27_week_26$phylogenetic_distance) #86.869
sd(c_scores_bcr_27_week_26$phylogenetic_distance) #36.176

#calculate nearest phylogenetic neighbor for each
bcr27week26_species <- colnames(bcr_27_week_26)
bcr27week26_min_phy_distances <- rep(NA, length(bcr27week26_species)) 

for(i in 1:length(bcr27week26_species)){
  bcr27week26_min_phy_distances[i] <- min(c_scores_bcr_27_week_26[c_scores_bcr_27_week_26$Sp1 == bcr27week26_species[i] |
                                                                    c_scores_bcr_27_week_26$Sp2 == bcr27week26_species[i], 
                                                                  "phylogenetic_distance"])
}

hist(bcr27week26_min_phy_distances)
mean(bcr27week26_min_phy_distances) #22.658
sd(bcr27week26_min_phy_distances) #15.878

hist(bcr27week26_min_phy_distances, freq=FALSE, col=rgb(0,0,1,1/2), 
     xlim=c(0, max(c_scores_bcr_27_week_26$phylogenetic_distance)), main="Phylogenetic Distances for NTI & NRI, BCR 27, Week 26",
     xlab="Phylogenetic Distance", breaks=20)
hist(c_scores_bcr_27_week_26$phylogenetic_distance, freq=FALSE, 
     breaks=20, col=rgb(1,0,0,1/2), add=TRUE)
abline(v=mean(c_scores_bcr_27_week_26$phylogenetic_distance), col="red", lwd=3)
abline(v=mean(bcr27week26_min_phy_distances), col="blue", lwd=3)



#BCR 27 week 1
hist(c_scores_bcr_27_week_1$phylogenetic_distance)
mean(c_scores_bcr_27_week_1$phylogenetic_distance) #77.379
sd(c_scores_bcr_27_week_1$phylogenetic_distance) #35.295

#calculate nearest phylogenetic neighbor for each
bcr27week1_species <- colnames(bcr_27_week_1)
bcr27week1_min_phy_distances <- rep(NA, length(bcr27week1_species)) 

for(i in 1:length(bcr27week1_species)){
  bcr27week1_min_phy_distances[i] <- min(c_scores_bcr_27_week_1[c_scores_bcr_27_week_1$Sp1 == bcr27week1_species[i] |
                                                                  c_scores_bcr_27_week_1$Sp2 == bcr27week1_species[i], 
                                                                "phylogenetic_distance"])
}

hist(bcr27week1_min_phy_distances)
mean(bcr27week1_min_phy_distances) #21.560
sd(bcr27week1_min_phy_distances) #17.896


hist(c_scores_bcr_27_week_1$phylogenetic_distance, freq=FALSE, 
     breaks=20, col=rgb(1,0,0,1/2), main="Phylogenetic Distances for NTI & NRI, BCR 27, Week 1",
     xlab="Phylogenetic Distance")
hist(bcr27week1_min_phy_distances, freq=FALSE, col=rgb(0,0,1,1/2), 
     xlim=c(0, max(c_scores_bcr_27_week_1$phylogenetic_distance)), breaks=20, add=TRUE)
abline(v=mean(c_scores_bcr_27_week_1$phylogenetic_distance), col="red", lwd=3)
abline(v=mean(bcr27week1_min_phy_distances), col="blue", lwd=3)


##### create a histogram with distributions of phylogenetic distance means across all BCRs
for(i in 5:37){
  assign(paste0("bcr_", i, "_week_26"), read.csv(paste0("output/bcr_community_data/26/community_26_bcr_", i, ".csv")))
} #takes about an hour

for(i in 5:37){
  assign(paste0("bcr_", i, "_week_1"), read.csv(paste0("output/bcr_community_data/1/community_1_bcr_", i, ".csv")))
} #takes about an hour

bcr_week_26_list <- list(bcr_5_week_26, bcr_6_week_26, bcr_7_week_26, bcr_8_week_26, bcr_9_week_26, bcr_10_week_26,
                         bcr_11_week_26, bcr_12_week_26, bcr_13_week_26, bcr_14_week_26, bcr_15_week_26, bcr_16_week_26,
                         bcr_17_week_26, bcr_18_week_26, bcr_19_week_26, bcr_20_week_26, bcr_21_week_26, bcr_22_week_26,
                         bcr_23_week_26, bcr_24_week_26, bcr_25_week_26, bcr_26_week_26, bcr_27_week_26, bcr_28_week_26,
                         bcr_29_week_26, bcr_30_week_26, bcr_31_week_26, bcr_32_week_26, bcr_33_week_26, bcr_34_week_26,
                         bcr_35_week_26, bcr_36_week_26, bcr_37_week_26)

#generate c scores for everyone at low frequencies so that we can generate phy distances
for(i in 5:37){
  df <- as.data.frame(bcr_week_26_list[i-4])
  df <- df[sample(nrow(df), nrow(df)/100, replace=FALSE),]
  ecospot <- ecospat.Cscore(df, 2, outpath="output/tables")
  assign(paste0("c_scores_bcr_", i, "_week_26"), read.delim("output/tables/Cscores.txt"))
}

#make a list of the c score matrices
c_scores_bcr_week_26_list <- list(c_scores_bcr_5_week_26, c_scores_bcr_6_week_26, c_scores_bcr_7_week_26, c_scores_bcr_8_week_26, c_scores_bcr_9_week_26, c_scores_bcr_10_week_26,
                                  c_scores_bcr_11_week_26, c_scores_bcr_12_week_26, c_scores_bcr_13_week_26, c_scores_bcr_14_week_26, c_scores_bcr_15_week_26, c_scores_bcr_16_week_26,
                                  c_scores_bcr_17_week_26, c_scores_bcr_18_week_26, c_scores_bcr_19_week_26, c_scores_bcr_20_week_26, c_scores_bcr_21_week_26, c_scores_bcr_22_week_26,
                                  c_scores_bcr_23_week_26, c_scores_bcr_24_week_26, c_scores_bcr_25_week_26, c_scores_bcr_26_week_26, c_scores_bcr_27_week_26, c_scores_bcr_28_week_26,
                                  c_scores_bcr_29_week_26, c_scores_bcr_30_week_26, c_scores_bcr_31_week_26, c_scores_bcr_32_week_26, c_scores_bcr_33_week_26, c_scores_bcr_34_week_26,
                                  c_scores_bcr_35_week_26, c_scores_bcr_36_week_26, c_scores_bcr_37_week_26)

#match phylogenies for each area

# Create a Phylogenetic Distance Matrix

#generate phylogeny matches
for(i in 5:37){
  assign(paste0("passerine.picCleanTree.bcr", i, "week26"), match.phylo.comm(phy = passerine.tree,
                                                                             comm = as.data.frame(bcr_week_26_list[i-4]))$phy)
  assign(paste0("passerine.picCleanComm.bcr", i, "week26"), match.phylo.comm(phy = passerine.tree,
                                                                             comm = as.data.frame(bcr_week_26_list[i-4]))$comm)
}

#make a list of those matches
passerine.picCleanTree.bcrweek_26_list <- list(passerine.picCleanTree.bcr5week26, passerine.picCleanTree.bcr6week26, passerine.picCleanTree.bcr7week26, passerine.picCleanTree.bcr8week26, passerine.picCleanTree.bcr9week26, passerine.picCleanTree.bcr10week26,
                                               passerine.picCleanTree.bcr11week26, passerine.picCleanTree.bcr12week26, passerine.picCleanTree.bcr13week26, passerine.picCleanTree.bcr14week26, passerine.picCleanTree.bcr15week26, passerine.picCleanTree.bcr16week26,
                                               passerine.picCleanTree.bcr17week26, passerine.picCleanTree.bcr18week26, passerine.picCleanTree.bcr19week26, passerine.picCleanTree.bcr20week26, passerine.picCleanTree.bcr21week26, passerine.picCleanTree.bcr22week26,
                                               passerine.picCleanTree.bcr23week26, passerine.picCleanTree.bcr24week26, passerine.picCleanTree.bcr25week26, passerine.picCleanTree.bcr26week26, passerine.picCleanTree.bcr27week26, passerine.picCleanTree.bcr28week26,
                                               passerine.picCleanTree.bcr29week26, passerine.picCleanTree.bcr30week26, passerine.picCleanTree.bcr31week26, passerine.picCleanTree.bcr32week26, passerine.picCleanTree.bcr33week26, passerine.picCleanTree.bcr34week26,
                                               passerine.picCleanTree.bcr35week26, passerine.picCleanTree.bcr36week26, passerine.picCleanTree.bcr37week26)

# create phylogenetic distance matrices and make a list of those matrices
for(i in 5:37){
  assign(paste0("passerine.cophenDist.bcr", i, "week26"), cophenetic.phylo(passerine.picCleanTree.bcrweek_26_list[[i-4]]))
}

passerine.cophenDist.week26_list <- list(passerine.cophenDist.bcr5week26, passerine.cophenDist.bcr6week26, passerine.cophenDist.bcr7week26, passerine.cophenDist.bcr8week26, passerine.cophenDist.bcr9week26, passerine.cophenDist.bcr10week26,
                                         passerine.cophenDist.bcr11week26, passerine.cophenDist.bcr12week26, passerine.cophenDist.bcr13week26, passerine.cophenDist.bcr14week26, passerine.cophenDist.bcr15week26, passerine.cophenDist.bcr16week26,
                                         passerine.cophenDist.bcr17week26, passerine.cophenDist.bcr18week26, passerine.cophenDist.bcr19week26, passerine.cophenDist.bcr20week26, passerine.cophenDist.bcr21week26, passerine.cophenDist.bcr22week26,
                                         passerine.cophenDist.bcr23week26, passerine.cophenDist.bcr24week26, passerine.cophenDist.bcr25week26, passerine.cophenDist.bcr26week26, passerine.cophenDist.bcr27week26, passerine.cophenDist.bcr28week26,
                                         passerine.cophenDist.bcr29week26, passerine.cophenDist.bcr30week26, passerine.cophenDist.bcr31week26, passerine.cophenDist.bcr32week26, passerine.cophenDist.bcr33week26, passerine.cophenDist.bcr34week26,
                                         passerine.cophenDist.bcr35week26, passerine.cophenDist.bcr36week26, passerine.cophenDist.bcr37week26)

#add phylogenetic distances to each of the c score matrices
for(i in 5:37){
  c_scores_bcr_week_26_list[[i-4]]$phylogenetic_distance <- rep(NA, nrow(c_scores_bcr_week_26_list[[i-4]]))
  for(j in 1:nrow(c_scores_bcr_week_26_list[[i-4]])){
    c_scores_bcr_week_26_list[[i-4]][j, "phylogenetic_distance"] <- passerine.cophenDist.week26_list[[i-4]][c_scores_bcr_week_26_list[[i-4]][j,"Sp1"],
                                                                                                            c_scores_bcr_week_26_list[[i-4]][j,"Sp2"]]
  }
}

#make a vector of mean NRI value within each BCR
week_26_mean_NRI_values <- rep(NA, length(c_scores_bcr_week_26_list))
for(i in 1:length(c_scores_bcr_week_26_list)){
  week_26_mean_NRI_values[i] <- mean(c_scores_bcr_week_26_list[[i]]$phylogenetic_distance)
}

#calculate nearest neighbors for each matrix
#generate species lists
for(i in 5:37){
  assign(paste0("bcr", i, "week26_species"), colnames(bcr_week_26_list[[i-4]]))
}

#generate list of species lists
week26_species_lists <- list(bcr5week26_species, bcr6week26_species, bcr7week26_species, bcr8week26_species, bcr9week26_species, bcr10week26_species,
                             bcr11week26_species, bcr12week26_species, bcr13week26_species, bcr14week26_species, bcr15week26_species, bcr16week26_species,
                             bcr17week26_species, bcr18week26_species, bcr19week26_species, bcr20week26_species, bcr21week26_species, bcr22week26_species,
                             bcr23week26_species, bcr24week26_species, bcr25week26_species, bcr26week26_species, bcr27week26_species, bcr28week26_species,
                             bcr29week26_species, bcr30week26_species, bcr31week26_species, bcr32week26_species, bcr33week26_species, bcr34week26_species,
                             bcr35week26_species, bcr36week26_species, bcr37week26_species)

#generate phy_distances matrices
for(i in 5:37){
  assign(paste0("bcr", i, "week26_min_phy_distances"), rep(NA, length(week26_species_lists[[i-4]])))
}

#generate list of phy_distance matrices
week26_min_phy_distances_lists <- list(bcr5week26_min_phy_distances, bcr6week26_min_phy_distances, bcr7week26_min_phy_distances, bcr8week26_min_phy_distances, bcr9week26_min_phy_distances, bcr10week26_min_phy_distances,
                                       bcr11week26_min_phy_distances, bcr12week26_min_phy_distances, bcr13week26_min_phy_distances, bcr14week26_min_phy_distances, bcr15week26_min_phy_distances, bcr16week26_min_phy_distances,
                                       bcr17week26_min_phy_distances, bcr18week26_min_phy_distances, bcr19week26_min_phy_distances, bcr20week26_min_phy_distances, bcr21week26_min_phy_distances, bcr22week26_min_phy_distances,
                                       bcr23week26_min_phy_distances, bcr24week26_min_phy_distances, bcr25week26_min_phy_distances, bcr26week26_min_phy_distances, bcr27week26_min_phy_distances, bcr28week26_min_phy_distances,
                                       bcr29week26_min_phy_distances, bcr30week26_min_phy_distances, bcr31week26_min_phy_distances, bcr32week26_min_phy_distances, bcr33week26_min_phy_distances, bcr34week26_min_phy_distances,
                                       bcr35week26_min_phy_distances, bcr36week26_min_phy_distances, bcr37week26_min_phy_distances)

#populate min phy distance matrices
for(i in 5:37){
  for(j in 1:length(week26_min_phy_distances_lists[[i-4]])){
    week26_min_phy_distances_lists[[i-4]][j] <- min(c_scores_bcr_week_26_list[[i-4]][c_scores_bcr_week_26_list[[i-4]]$Sp1 == week26_species_lists[[i-4]][j] |
                                                                                       c_scores_bcr_week_26_list[[i-4]]$Sp2 == week26_species_lists[[i-4]][j],
                                                                                     "phylogenetic_distance"])
  }
}

#make a vector of mean NTI value within each BCR
week_26_mean_NTI_values <- rep(NA, length(c_scores_bcr_week_26_list))
for(i in 1:length(c_scores_bcr_week_26_list)){
  week_26_mean_NTI_values[i] <- mean(week26_min_phy_distances_lists[[i]])
}

#initial histogram of mean NRI and NTI values within each BCR
hist(week_26_mean_NRI_values, freq=FALSE, col=rgb(0,0,1,1/2), 
     xlim=c(0, max(week_26_mean_NRI_values)), ylim=c(0,.25), main="BCR Mean Phylogenetic Distances for NTI & NRI Week 26",
     xlab="Phylogenetic Distance", breaks=10)
hist(week_26_mean_NTI_values, freq=FALSE, 
     breaks=5, col=rgb(1,0,0,1/2), add=TRUE)
abline(v=mean(week_26_mean_NRI_values), col="blue", lwd=2)
abline(v=mean(week_26_mean_NTI_values), col="red", lwd=2)


#### replicate process for week 1
bcr_week_1_list <- list(bcr_5_week_1, bcr_6_week_1, bcr_7_week_1, bcr_8_week_1, bcr_9_week_1, bcr_10_week_1,
                        bcr_11_week_1, bcr_12_week_1, bcr_13_week_1, bcr_14_week_1, bcr_15_week_1, bcr_16_week_1,
                        bcr_17_week_1, bcr_18_week_1, bcr_19_week_1, bcr_20_week_1, bcr_21_week_1, bcr_22_week_1,
                        bcr_23_week_1, bcr_24_week_1, bcr_25_week_1, bcr_26_week_1, bcr_27_week_1, bcr_28_week_1,
                        bcr_29_week_1, bcr_30_week_1, bcr_31_week_1, bcr_32_week_1, bcr_33_week_1, bcr_34_week_1,
                        bcr_35_week_1, bcr_36_week_1, bcr_37_week_1)

#generate c scores for everyone at low frequencies so that we can generate phy distances
for(i in 5:37){
  df <- as.data.frame(bcr_week_1_list[i-4])
  df <- df[sample(nrow(df), nrow(df)/100, replace=FALSE),]
  ecospot <- ecospat.Cscore(df, 2, outpath="output/tables")
  assign(paste0("c_scores_bcr_", i, "_week_1"), read.delim("output/tables/Cscores.txt"))
}

#make a list of the c score matrices
c_scores_bcr_week_1_list <- list(c_scores_bcr_5_week_1, c_scores_bcr_6_week_1, c_scores_bcr_7_week_1, c_scores_bcr_8_week_1, c_scores_bcr_9_week_1, c_scores_bcr_10_week_1,
                                 c_scores_bcr_11_week_1, c_scores_bcr_12_week_1, c_scores_bcr_13_week_1, c_scores_bcr_14_week_1, c_scores_bcr_15_week_1, c_scores_bcr_16_week_1,
                                 c_scores_bcr_17_week_1, c_scores_bcr_18_week_1, c_scores_bcr_19_week_1, c_scores_bcr_20_week_1, c_scores_bcr_21_week_1, c_scores_bcr_22_week_1,
                                 c_scores_bcr_23_week_1, c_scores_bcr_24_week_1, c_scores_bcr_25_week_1, c_scores_bcr_26_week_1, c_scores_bcr_27_week_1, c_scores_bcr_28_week_1,
                                 c_scores_bcr_29_week_1, c_scores_bcr_30_week_1, c_scores_bcr_31_week_1, c_scores_bcr_32_week_1, c_scores_bcr_33_week_1, c_scores_bcr_34_week_1,
                                 c_scores_bcr_35_week_1, c_scores_bcr_36_week_1, c_scores_bcr_37_week_1)

#match phylogenies for each area

# Create a Phylogenetic Distance Matrix

#generate phylogeny matches
for(i in 5:37){
  assign(paste0("passerine.picCleanTree.bcr", i, "week1"), match.phylo.comm(phy = passerine.tree,
                                                                            comm = as.data.frame(bcr_week_1_list[i-4]))$phy)
  assign(paste0("passerine.picCleanComm.bcr", i, "week1"), match.phylo.comm(phy = passerine.tree,
                                                                            comm = as.data.frame(bcr_week_1_list[i-4]))$comm)
}

#make a list of those matches
passerine.picCleanTree.bcrweek_1_list <- list(passerine.picCleanTree.bcr5week1, passerine.picCleanTree.bcr6week1, passerine.picCleanTree.bcr7week1, passerine.picCleanTree.bcr8week1, passerine.picCleanTree.bcr9week1, passerine.picCleanTree.bcr10week1,
                                              passerine.picCleanTree.bcr11week1, passerine.picCleanTree.bcr12week1, passerine.picCleanTree.bcr13week1, passerine.picCleanTree.bcr14week1, passerine.picCleanTree.bcr15week1, passerine.picCleanTree.bcr16week1,
                                              passerine.picCleanTree.bcr17week1, passerine.picCleanTree.bcr18week1, passerine.picCleanTree.bcr19week1, passerine.picCleanTree.bcr20week1, passerine.picCleanTree.bcr21week1, passerine.picCleanTree.bcr22week1,
                                              passerine.picCleanTree.bcr23week1, passerine.picCleanTree.bcr24week1, passerine.picCleanTree.bcr25week1, passerine.picCleanTree.bcr26week1, passerine.picCleanTree.bcr27week1, passerine.picCleanTree.bcr28week1,
                                              passerine.picCleanTree.bcr29week1, passerine.picCleanTree.bcr30week1, passerine.picCleanTree.bcr31week1, passerine.picCleanTree.bcr32week1, passerine.picCleanTree.bcr33week1, passerine.picCleanTree.bcr34week1,
                                              passerine.picCleanTree.bcr35week1, passerine.picCleanTree.bcr36week1, passerine.picCleanTree.bcr37week1)

# create phylogenetic distance matrices and make a list of those matrices
for(i in 5:37){
  assign(paste0("passerine.cophenDist.bcr", i, "week1"), cophenetic.phylo(passerine.picCleanTree.bcrweek_1_list[[i-4]]))
}

passerine.cophenDist.week1_list <- list(passerine.cophenDist.bcr5week1, passerine.cophenDist.bcr6week1, passerine.cophenDist.bcr7week1, passerine.cophenDist.bcr8week1, passerine.cophenDist.bcr9week1, passerine.cophenDist.bcr10week1,
                                        passerine.cophenDist.bcr11week1, passerine.cophenDist.bcr12week1, passerine.cophenDist.bcr13week1, passerine.cophenDist.bcr14week1, passerine.cophenDist.bcr15week1, passerine.cophenDist.bcr16week1,
                                        passerine.cophenDist.bcr17week1, passerine.cophenDist.bcr18week1, passerine.cophenDist.bcr19week1, passerine.cophenDist.bcr20week1, passerine.cophenDist.bcr21week1, passerine.cophenDist.bcr22week1,
                                        passerine.cophenDist.bcr23week1, passerine.cophenDist.bcr24week1, passerine.cophenDist.bcr25week1, passerine.cophenDist.bcr26week1, passerine.cophenDist.bcr27week1, passerine.cophenDist.bcr28week1,
                                        passerine.cophenDist.bcr29week1, passerine.cophenDist.bcr30week1, passerine.cophenDist.bcr31week1, passerine.cophenDist.bcr32week1, passerine.cophenDist.bcr33week1, passerine.cophenDist.bcr34week1,
                                        passerine.cophenDist.bcr35week1, passerine.cophenDist.bcr36week1, passerine.cophenDist.bcr37week1)

#add phylogenetic distances to each of the c score matrices
for(i in 5:37){
  c_scores_bcr_week_1_list[[i-4]]$phylogenetic_distance <- rep(NA, nrow(c_scores_bcr_week_1_list[[i-4]]))
  for(j in 1:nrow(c_scores_bcr_week_1_list[[i-4]])){
    c_scores_bcr_week_1_list[[i-4]][j, "phylogenetic_distance"] <- passerine.cophenDist.week1_list[[i-4]][c_scores_bcr_week_1_list[[i-4]][j,"Sp1"],
                                                                                                          c_scores_bcr_week_1_list[[i-4]][j,"Sp2"]]
  }
}

#make a vector of mean NRI value within each BCR
week_1_mean_NRI_values <- rep(NA, length(c_scores_bcr_week_1_list))
for(i in 1:length(c_scores_bcr_week_1_list)){
  week_1_mean_NRI_values[i] <- mean(c_scores_bcr_week_1_list[[i]]$phylogenetic_distance)
}

#calculate nearest neighbors for each matrix
#generate species lists
for(i in 5:37){
  assign(paste0("bcr", i, "week1_species"), colnames(bcr_week_1_list[[i-4]]))
}

#generate list of species lists
week1_species_lists <- list(bcr5week1_species, bcr6week1_species, bcr7week1_species, bcr8week1_species, bcr9week1_species, bcr10week1_species,
                            bcr11week1_species, bcr12week1_species, bcr13week1_species, bcr14week1_species, bcr15week1_species, bcr16week1_species,
                            bcr17week1_species, bcr18week1_species, bcr19week1_species, bcr20week1_species, bcr21week1_species, bcr22week1_species,
                            bcr23week1_species, bcr24week1_species, bcr25week1_species, bcr26week1_species, bcr27week1_species, bcr28week1_species,
                            bcr29week1_species, bcr30week1_species, bcr31week1_species, bcr32week1_species, bcr33week1_species, bcr34week1_species,
                            bcr35week1_species, bcr36week1_species, bcr37week1_species)

#generate phy_distances matrices
for(i in 5:37){
  assign(paste0("bcr", i, "week1_min_phy_distances"), rep(NA, length(week1_species_lists[[i-4]])))
}

#generate list of phy_distance matrices
week1_min_phy_distances_lists <- list(bcr5week1_min_phy_distances, bcr6week1_min_phy_distances, bcr7week1_min_phy_distances, bcr8week1_min_phy_distances, bcr9week1_min_phy_distances, bcr10week1_min_phy_distances,
                                      bcr11week1_min_phy_distances, bcr12week1_min_phy_distances, bcr13week1_min_phy_distances, bcr14week1_min_phy_distances, bcr15week1_min_phy_distances, bcr16week1_min_phy_distances,
                                      bcr17week1_min_phy_distances, bcr18week1_min_phy_distances, bcr19week1_min_phy_distances, bcr20week1_min_phy_distances, bcr21week1_min_phy_distances, bcr22week1_min_phy_distances,
                                      bcr23week1_min_phy_distances, bcr24week1_min_phy_distances, bcr25week1_min_phy_distances, bcr26week1_min_phy_distances, bcr27week1_min_phy_distances, bcr28week1_min_phy_distances,
                                      bcr29week1_min_phy_distances, bcr30week1_min_phy_distances, bcr31week1_min_phy_distances, bcr32week1_min_phy_distances, bcr33week1_min_phy_distances, bcr34week1_min_phy_distances,
                                      bcr35week1_min_phy_distances, bcr36week1_min_phy_distances, bcr37week1_min_phy_distances)

#populate min phy distance matrices
for(i in 5:37){
  for(j in 1:length(week1_min_phy_distances_lists[[i-4]])){
    week1_min_phy_distances_lists[[i-4]][j] <- min(c_scores_bcr_week_1_list[[i-4]][c_scores_bcr_week_1_list[[i-4]]$Sp1 == week1_species_lists[[i-4]][j] |
                                                                                     c_scores_bcr_week_1_list[[i-4]]$Sp2 == week1_species_lists[[i-4]][j],
                                                                                   "phylogenetic_distance"])
  }
}

#make a vector of mean NTI value within each BCR
week_1_mean_NTI_values <- rep(NA, length(c_scores_bcr_week_1_list))
for(i in 1:length(c_scores_bcr_week_1_list)){
  week_1_mean_NTI_values[i] <- mean(week1_min_phy_distances_lists[[i]])
}

#initial histogram of mean NRI and NTI values within each BCR
hist(week_1_mean_NRI_values, freq=FALSE, col=rgb(0,0,1,1/2), 
     xlim=c(0, max(week_1_mean_NRI_values)), ylim=c(0,.25), main="BCR Mean Phylogenetic Distances for NTI & NRI Week 1",
     xlab="Phylogenetic Distance", breaks=10)
hist(week_1_mean_NTI_values, freq=FALSE, 
     breaks=5, col=rgb(1,0,0,1/2), add=TRUE)
abline(v=mean(week_1_mean_NRI_values), col="blue", lwd=2)
abline(v=mean(week_1_mean_NTI_values), col="red", lwd=2)


#### create nicer density plots for manuscript
#week 26
mean_NRI_NTI_hist_26 <- data.frame(matrix(NA, nrow = 2*length(week_26_mean_NRI_values), ncol = 2))
colnames(mean_NRI_NTI_hist_26) <- c("value", "type")
mean_NRI_NTI_hist_26[c(1:length(week_26_mean_NRI_values)), "value"] <- week_26_mean_NRI_values
mean_NRI_NTI_hist_26[c((length(week_26_mean_NRI_values) + 1)):(nrow(mean_NRI_NTI_hist_26)), "value"] <- week_26_mean_NTI_values
mean_NRI_NTI_hist_26[c(1:length(week_26_mean_NRI_values)), "type"] <- "NRI"
mean_NRI_NTI_hist_26[c((length(week_26_mean_NRI_values) + 1)):(nrow(mean_NRI_NTI_hist_26)), "type"] <- "NTI"
mean_NRI_NTI_hist_26

week26density <- ggplot(mean_NRI_NTI_hist_26, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("Phylogenetic Distance") +
  ylab("Density") + 
  xlim(0,100) +
  guides(fill=guide_legend(title="Metric")) +
  ggtitle("Week 26") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25))

#week 1
mean_NRI_NTI_hist_1 <- data.frame(matrix(NA, nrow = 2*length(week_1_mean_NRI_values), ncol = 2))
colnames(mean_NRI_NTI_hist_1) <- c("value", "type")
mean_NRI_NTI_hist_1[c(1:length(week_1_mean_NRI_values)), "value"] <- week_1_mean_NRI_values
mean_NRI_NTI_hist_1[c((length(week_1_mean_NRI_values) + 1)):(nrow(mean_NRI_NTI_hist_1)), "value"] <- week_1_mean_NTI_values
mean_NRI_NTI_hist_1[c(1:length(week_1_mean_NRI_values)), "type"] <- "NRI"
mean_NRI_NTI_hist_1[c((length(week_1_mean_NRI_values) + 1)):(nrow(mean_NRI_NTI_hist_1)), "type"] <- "NTI"
mean_NRI_NTI_hist_1

week1density <- ggplot(mean_NRI_NTI_hist_1, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("Phylogenetic Distance") +
  ylab("Density") + 
  xlim(0,100) +
  guides(fill=guide_legend(title="Metric")) +
  ggtitle("Week 1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25))

#merge two dataframes
both_weeks_density <- dplyr::bind_rows(mean_NRI_NTI_hist_26, mean_NRI_NTI_hist_1, .id = "week")
both_weeks_density[both_weeks_density$week == "2","week"] <- "Week 1"
both_weeks_density[both_weeks_density$week == "1","week"] <- "Week 26"

library(ggpubr)
ggarrange(week1density, week26density, ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")

ggplot(both_weeks_density, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("Phylogenetic Distance") +
  ylab("Density") + 
  xlim(0,100) +
  guides(fill=guide_legend(title="Metric")) +
  ggtitle("Mean NRI vs. NTI Phylogenetic Distances") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~week) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25)) 

mean(both_weeks_density[both_weeks_density$week == "Week 26" & both_weeks_density$type == "NRI", "value"]) #week 26 NRI: 86.2577
mean(both_weeks_density[both_weeks_density$week == "Week 26" & both_weeks_density$type == "NTI", "value"]) #week 26 NTI: 19.6550; ratio = 4.386
mean(both_weeks_density[both_weeks_density$week == "Week 1" & both_weeks_density$type == "NRI", "value"]) #week 1 NRI: 80.1144
mean(both_weeks_density[both_weeks_density$week == "Week 1" & both_weeks_density$type == "NTI", "value"]) #week 1 NTI: 22.6048; ratio = 3.544


##### some more experimentation with turnover
#(load all_combined_NRI_df and all_combined_NTI_df)
community_intersect_26_1 <- intersect(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26,"X"],
                                      all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1,"X"]) #447 intersecting

NTI_week_1 <- rep(NA, length(community_intersect_26_1))

all_combined_week_NTIs_df[all_combined_week_NTIs_df$X == community_intersect_26_1[1] &
                            all_combined_week_NTIs_df$.id == 1, "NTI"]

NTI_turnover_ratios <- data.frame(community_intersect_26_1, NTI_week_1)
NTI_turnover_ratios$NTI_week_26 <- rep(NA, length(community_intersect_26_1))
NTI_turnover_ratios$richness_week_1 <- rep(NA, length(community_intersect_26_1))
NTI_turnover_ratios$richness_week_26 <- rep(NA, length(community_intersect_26_1))
NTI_turnover_ratios$NTI_turnover <- rep(NA, length(community_intersect_26_1))
NTI_turnover_ratios$richness_turnover <- rep(NA, length(community_intersect_26_1))
NTI_turnover_ratios$BCR <- rep(NA, length(community_intersect_26_1))
NTI_turnover_ratios$mean_latitude <- rep(NA, length(community_intersect_26_1))
NTI_turnover_ratios$latitude <- rep(NA, length(community_intersect_26_1))
NTI_turnover_ratios$longitude <- rep(NA, length(community_intersect_26_1))
NTI_turnover_ratios$mig_prop_week_1 <- rep(NA, length(community_intersect_26_1))
NTI_turnover_ratios$mig_prop_week_26 <- rep(NA, length(community_intersect_26_1))
NTI_turnover_ratios$mig_prop_turnover <- rep(NA, length(community_intersect_26_1))



mean_latitudes <- c(50.6826, 58.50326, 57.83681, 54.43946, 43.24103, 46.19079,
                    47.40598, 47.95409, 43.65609, 44.72067, 38.88149, 38.52856,
                    44.91282, 38.69052, 38.43332, 30.43795, 34.65978, 40.8388,
                    42.93611, 37.05901, 33.09216, 35.14722, 33.5757, 37.8537,
                    37.49024, 40.89619, 27.68041, 34.09171, 32.76825, 25.70883,
                    26.28559, 26.81768, 28.57299)


for(i in 1:nrow(NTI_turnover_ratios)){
  NTI_turnover_ratios[i, "NTI_week_1"] <- all_combined_week_NTIs_df[all_combined_week_NTIs_df$X == community_intersect_26_1[i] &
                                                                      all_combined_week_NTIs_df$.id == 1, "NTI"]
  NTI_turnover_ratios[i, "NTI_week_26"] <- all_combined_week_NTIs_df[all_combined_week_NTIs_df$X == community_intersect_26_1[i] &
                                                                       all_combined_week_NTIs_df$.id == 26, "NTI"]
  NTI_turnover_ratios[i, "richness_week_1"] <- all_combined_week_NTIs_df[all_combined_week_NTIs_df$X == community_intersect_26_1[i] &
                                                                           all_combined_week_NTIs_df$.id == 1, "ntaxa"]
  NTI_turnover_ratios[i, "richness_week_26"] <- all_combined_week_NTIs_df[all_combined_week_NTIs_df$X == community_intersect_26_1[i] &
                                                                            all_combined_week_NTIs_df$.id == 26, "ntaxa"]
  NTI_turnover_ratios[i, "BCR"] <- all_combined_week_NTIs_df[all_combined_week_NTIs_df$X == community_intersect_26_1[i], "BCR"][1]
  NTI_turnover_ratios[i, "latitude"] <- all_combined_week_NTIs_df[all_combined_week_NTIs_df$X == community_intersect_26_1[i], "latitude"][1]
  NTI_turnover_ratios[i, "longitude"] <- all_combined_week_NTIs_df[all_combined_week_NTIs_df$X == community_intersect_26_1[i], "longitude"][1]
  NTI_turnover_ratios[i, "mig_prop_week_1"] <- all_combined_week_NTIs_df[all_combined_week_NTIs_df$X == community_intersect_26_1[i] &
                                                                           all_combined_week_NTIs_df$.id == 1, "proportion"]
  NTI_turnover_ratios[i, "mig_prop_week_26"] <- all_combined_week_NTIs_df[all_combined_week_NTIs_df$X == community_intersect_26_1[i] &
                                                                           all_combined_week_NTIs_df$.id == 26, "proportion"]
}

for(i in 1:nrow(NTI_turnover_ratios)){
  NTI_turnover_ratios[i, "NTI_turnover"] <- NTI_turnover_ratios[i, "NTI_week_26"] - NTI_turnover_ratios[i, "NTI_week_1"]
  NTI_turnover_ratios[i, "richness_turnover"] <- NTI_turnover_ratios[i, "richness_week_26"] - NTI_turnover_ratios[i, "richness_week_1"]
  NTI_turnover_ratios[i, "mig_prop_turnover"] <- NTI_turnover_ratios[i, "mig_prop_week_26"] - NTI_turnover_ratios[i, "mig_prop_week_1"]
  NTI_turnover_ratios[i, "mean_latitude"] <- mean_latitudes[NTI_turnover_ratios[i, "BCR"]-4]
}  

plot(NTI_turnover_ratios$NTI_turnover ~ NTI_turnover_ratios$richness_turnover) #hella significant


### NRI turnover experiment

community_intersect_26_1_NRI <- intersect(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 26,"X"],
                                          all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 1,"X"]) #447 intersecting

NRI_week_1 <- rep(NA, length(community_intersect_26_1_NRI))

all_combined_week_NRIs_df[all_combined_week_NRIs_df$X == community_intersect_26_1_NRI[1] &
                            all_combined_week_NRIs_df$.id == 1, "NRI"]

NRI_turnover_ratios <- data.frame(community_intersect_26_1_NRI, NRI_week_1)
NRI_turnover_ratios$NRI_week_26 <- rep(NA, length(community_intersect_26_1_NRI))
NRI_turnover_ratios$richness_week_1 <- rep(NA, length(community_intersect_26_1_NRI))
NRI_turnover_ratios$richness_week_26 <- rep(NA, length(community_intersect_26_1_NRI))
NRI_turnover_ratios$NRI_turnover <- rep(NA, length(community_intersect_26_1_NRI))
NRI_turnover_ratios$richness_turnover <- rep(NA, length(community_intersect_26_1_NRI))
NRI_turnover_ratios$BCR <- rep(NA, length(community_intersect_26_1))
NRI_turnover_ratios$mean_latitude <- rep(NA, length(community_intersect_26_1))
NRI_turnover_ratios$latitude <- rep(NA, length(community_intersect_26_1))
NRI_turnover_ratios$longitude <- rep(NA, length(community_intersect_26_1))
NRI_turnover_ratios$mig_prop_week_1 <- rep(NA, length(community_intersect_26_1))
NRI_turnover_ratios$mig_prop_week_26 <- rep(NA, length(community_intersect_26_1))
NRI_turnover_ratios$mig_prop_turnover <- rep(NA, length(community_intersect_26_1))

for(i in 1:nrow(NRI_turnover_ratios)){
  NRI_turnover_ratios[i, "NRI_week_1"] <- all_combined_week_NRIs_df[all_combined_week_NRIs_df$X == community_intersect_26_1_NRI[i] &
                                                                      all_combined_week_NRIs_df$.id == 1, "NRI"]
  NRI_turnover_ratios[i, "NRI_week_26"] <- all_combined_week_NRIs_df[all_combined_week_NRIs_df$X == community_intersect_26_1_NRI[i] &
                                                                       all_combined_week_NRIs_df$.id == 26, "NRI"]
  NRI_turnover_ratios[i, "richness_week_1"] <- all_combined_week_NRIs_df[all_combined_week_NRIs_df$X == community_intersect_26_1_NRI[i] &
                                                                           all_combined_week_NRIs_df$.id == 1, "ntaxa"]
  NRI_turnover_ratios[i, "richness_week_26"] <- all_combined_week_NRIs_df[all_combined_week_NRIs_df$X == community_intersect_26_1_NRI[i] &
                                                                            all_combined_week_NRIs_df$.id == 26, "ntaxa"]
  NRI_turnover_ratios[i, "BCR"] <- all_combined_week_NRIs_df[all_combined_week_NRIs_df$X == community_intersect_26_1[i], "BCR"][1]
  NRI_turnover_ratios[i, "latitude"] <- all_combined_week_NRIs_df[all_combined_week_NRIs_df$X == community_intersect_26_1[i], "latitude"][1]
  NRI_turnover_ratios[i, "longitude"] <- all_combined_week_NRIs_df[all_combined_week_NRIs_df$X == community_intersect_26_1[i], "longitude"][1]
  NRI_turnover_ratios[i, "mig_prop_week_1"] <- all_combined_week_NRIs_df[all_combined_week_NRIs_df$X == community_intersect_26_1[i] &
                                                                           all_combined_week_NRIs_df$.id == 1, "proportion"]
  NRI_turnover_ratios[i, "mig_prop_week_26"] <- all_combined_week_NRIs_df[all_combined_week_NRIs_df$X == community_intersect_26_1[i] &
                                                                            all_combined_week_NRIs_df$.id == 26, "proportion"]
}

for(i in 1:nrow(NRI_turnover_ratios)){
  NRI_turnover_ratios[i, "NRI_turnover"] <- NRI_turnover_ratios[i, "NRI_week_26"] - NRI_turnover_ratios[i, "NRI_week_1"]
  NRI_turnover_ratios[i, "richness_turnover"] <- NRI_turnover_ratios[i, "richness_week_26"] - NRI_turnover_ratios[i, "richness_week_1"]
  NRI_turnover_ratios[i, "mig_prop_turnover"] <- NRI_turnover_ratios[i, "mig_prop_week_26"] - NRI_turnover_ratios[i, "mig_prop_week_1"]
  NRI_turnover_ratios[i, "mean_latitude"] <- mean_latitudes[NRI_turnover_ratios[i, "BCR"]-4]
}  

plot(NRI_turnover_ratios$NRI_turnover ~ NRI_turnover_ratios$richness_turnover) ## also hella significant!


library(ggExtra)
#### first try to plot them nicer
#plot NTI turnover vs. richness turnover
NTI_turnover_plot <- ggplot(mapping = aes(x = NTI_turnover_ratios$richness_turnover, y = NTI_turnover_ratios$NTI_turnover)) + 
  geom_point(aes(color = NTI_turnover_ratios$latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("NTI vs. Richness Differences (June - January)") + 
  labs(x = "Richness Difference", y = "NTI Difference") +
  geom_smooth(method = "lm")
ggMarginal(NTI_turnover_plot, fill = "green")

NTI_turnover_gls <- gls(NTI_turnover ~ richness_turnover, data = NTI_turnover_ratios, corSpher(form = ~latitude+longitude, nugget = TRUE))
summary(NTI_turnover_gls) #p = 0.0019

#plot NRI turnover vs. richness turnover
ggplot(mapping = aes(x = NRI_turnover_ratios$richness_turnover, y = NRI_turnover_ratios$NRI_turnover)) + 
  geom_point(aes(color = NRI_turnover_ratios$latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("NRI vs. Richness Differences (June - January)") + 
  labs(x = "Richness Difference", y = "NRI Difference") +
  geom_smooth(method = "lm")

NRI_turnover_gls <- gls(NRI_turnover ~ richness_turnover, data = NRI_turnover_ratios, corSpher(form = ~latitude+longitude, nugget = TRUE))
summary(NRI_turnover_gls) #p = 0

#plot NTI turnover vs. mig prop turnover
NTI_mig_prop_turnover_plot <- ggplot(mapping = aes(x = NTI_turnover_ratios$mig_prop_turnover, y = NTI_turnover_ratios$NTI_turnover)) + 
  geom_point(aes(color = NTI_turnover_ratios$latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("NTI vs. Mig. Prop. Differences (June - January)") + 
  labs(x = "Mig. Prop. Difference", y = "NTI Difference") +
  geom_smooth(method = "lm")
ggMarginal(NTI_mig_prop_turnover_plot, fill = "green")


NTI_mig_prop_turnover_gls <- gls(NTI_turnover ~ mig_prop_turnover, data = NTI_turnover_ratios, corSpher(form = ~latitude+longitude, nugget = TRUE))
summary(NTI_mig_prop_turnover_gls) #p = 0


#plot NRI turnover vs. mig prop turnover
NRI_mig_prop_turnover_plot <- ggplot(mapping = aes(x = NRI_turnover_ratios$mig_prop_turnover, y = NRI_turnover_ratios$NRI_turnover)) + 
  geom_point(aes(color = NRI_turnover_ratios$latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("NRI vs. Mig. Prop. Differences (June - January)") + 
  labs(x = "Mig. Prop. Difference", y = "NRI Difference") +
  geom_smooth(method = "lm")
ggMarginal(NRI_mig_prop_turnover_plot, fill = "green")


NRI_mig_prop_turnover_gls <- gls(NRI_turnover ~ mig_prop_turnover, data = NRI_turnover_ratios, corSpher(form = ~latitude+longitude, nugget = TRUE))
summary(NRI_mig_prop_turnover_gls) #p = 0


### replicate with continental
### reload 




#### now model with statistical corrections
#library(ecodados)
#library(visdat)
library(tidyverse)
library(lattice)
library(RVAideMemoire)
library(DHARMa)
library(performance)
#library(MuMIn)
#library(piecewiseSEM)
library(MASS)
library(ggExtra)
library(Rmisc)
library(emmeans) 
library(sjPlot)
library(bbmle)
library(glmmTMB)
library(ordinal)
library(car)
library(ecolottery)
library(naniar)
library(vcd)
library(generalhoslem)


