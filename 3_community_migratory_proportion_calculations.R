# 3. Community migratory proportion calculations
# This script defines a function that counts the number of migrant species in each community.
# The proportion of migratory species in each community is then calculated.
# These proportions are stored with other community phylogenetic metrics.
# Script gives an example for week 1, but can be replicated for other 51 weeks of the year with minimal substitutions (indicated below) 

#load required libraries
required.libraries <- c("ebirdst","terra","sf","dplyr","tidyr","rnaturalearth",
                        "geodata","ggplot2","fields","stringr","TreeTools",
                        "viridis","rgdal","raster","devtools","data.table",
                        "ape","picante","pez","phytools","vegan","adephylo",
                        "phylobase","geiger","mvMORPH","OUwie","hisse",
                        "BAMMtools","phylosignal","Biostrings","kableExtra",
                        "betapart","gridExtra","reshape2")

lapply(required.libraries, require, character.only = TRUE)

project <- terra::project 

# load the mapping parameters to generate necessary CRS
path <- get_species_path("example_data")
fac_parameters <- load_fac_map_parameters(path)
crs <- fac_parameters$custom_projection

#read in CSV defining species codes and scientific names
species_codes <- read.csv("species_codes.csv")

# load list of species names to be included in dataset
species_names <- species_codes$ebirdst_code

#create vector of scientific names of the above species for interfacing with phylogenetic tree
species_passerines <- species_codes$jetz_tree_latin_name

#load in the phylogenetic tree
tree.passerine <- read.nexus("phylogeny_data/jetz_passerine_tree/jetz_consensus_tree.tre")
#define outgroup
tree.passerine <- root(tree.passerine, outgroup = c(263:297))
#visualize tree to make sure it is loaded in properly
plot(tree.passerine)

#drop tips to ensure that species in .tre file and list of species names match
passerine.tree <- drop.tip(tree.passerine,tree.passerine$tip.label[-match(species_passerines, tree.passerine$tip.label)])

#load example raster
ovenbi1_abd_seasonal_mean <- load_raster("data/ebirdst/2021/ovenbi1", product = "abundance", 
                                         period = "seasonal", metric = "mean", 
                                         resolution = "hr")
#use example raster to crop to US boundaries
us <- ne_countries(country = "united states of america", type = "countries", returnclass = "sf") %>% 
  st_transform(st_crs(ovenbi1_abd_seasonal_mean))

##### ASSIGN MIGRATORY PROPORTION VALUES TO EACH COMMUNITY #####

#use week 1 as an example; week values and community files will be replaced when calculating for other weeks

#load first taxonomic species at minimum community size (established in raster_download.R)
wewpew_1_agg_df <- read.csv("output/species_data/1/wewpew_1_agg.csv", row.names = 1)
#read existing community table (generated in raster_download.R)
community1 <- data.table::fread("output/community_data/community_1.csv")
#remove "V1" columns if they have been generated in loading in community data
community1 <- data.frame(community1[,-c("V1")])
#assign rownames to match the individually identified communities
rownames(community1) <- rownames(wewpew_1_agg_df)

#read in CSV defining migratory status for each species
migratory_status <- read.csv("migratory_status.csv")

#create empty dataframe to store values for number of migrants, number of species, and proportion of migrants in community
migratory_proportion_df_week_1 <- data.frame(community = rep(NA, nrow(community1)),
                                             nmigrant = rep(NA, nrow(community1)), #number of migrant species
                                             nspecies = rep(NA, nrow(community1)), #number of total species
                                             proportion = rep(NA, nrow(community1))) #proportion of species that are migratory

#populate community individual identifiers
migratory_proportion_df_week_1[,"community"] <- rownames(community1)

#populate species richness into nspecies column
migratory_proportion_df_week_1[,"nspecies"] <- apply(community1, 1, function(i) sum(i > 0))

#calculate number of migrants in each community
#create function for number of migrants in a community
number_of_migrants <- function(i){
  row_comm_present_species <- i
  test_migratory_status <- rep(NA, ncol(row_comm_present_species))
  for(i in 1:ncol(row_comm_present_species)){
    test_migratory_status[i] <- migratory_status[migratory_status$jetz_tree_latin_name == colnames(row_comm_present_species)[i],"migratory_status"]
  }
  return(length(test_migratory_status[test_migratory_status == "migratory"]))
}

#fill into dataset (warning: can take multiple hours)
for(i in 1:nrow(migratory_proportion_df_week_1)){
  migratory_proportion_df_week_1[i,"nmigrant"] <- number_of_migrants(community1[i,] %>% select_if(~any(. > 0)))
}

#calculate migratory proportion
migratory_proportion_df_week_1$proportion <- migratory_proportion_df_week_1$nmigrant/migratory_proportion_df_week_1$nspecies

#export file of migratory proportion values for each week 1 community
write.csv(migratory_proportion_df_week_1, "output/community_data/migratory_proportion_week_1.csv")

#reload in migratory proportion dataframes if already created
migratory_proportion_df_week_1 <- read.csv("output/community_data/migratory_proportion_week_1.csv")
migratory_proportion_df_week_1 <- data.frame(migratory_proportion_df_week_1[, c("community", "nmigrant", "nspecies", "proportion")],
                                             row.names = migratory_proportion_df_week_1$community) 

##### PROCESS MIGRATORY PROPORTION IN BCR-DIVIDED DATAFRAMES #####

#read in shapefiles of Bird Conservation Regions (BCRs)
allbcrs <- read_sf("data/BCR_Terrestrial/BCR_Terrestrial_master.shp")
#assign row names to BCRs
wewpew_1_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 1), us), crs, method = "near"), fact=2)
for(i in 1:37){
  assign(paste0("bcr",i), allbcrs[allbcrs$BCR == i,] %>%
           st_transform(crs = crs) %>%
           st_geometry())
  assign(paste0("bcr", i, "spatial"), as_Spatial(get(paste0("bcr", i))))
  assign(paste0("wewpew_1_bcr_", i), rast(mask(raster(wewpew_1_agg), get(paste0("bcr", i, "spatial")))))
  assign(paste0("wewpew_1_bcr_", i, "_df"), as.data.frame(get(paste0("wewpew_1_bcr_", i))))
  assign(paste0("bcr", i, "rows"), row.names(get(paste0("wewpew_1_bcr_", i, "_df"))))
}

#extract by BCR
migratory_proportion_df_week_1 <- data.frame(migratory_proportion_df_week_1[, c("community", "nmigrant", "nspecies", "proportion")],
                                             row.names = migratory_proportion_df_week_1$community) 

#subset migratory proportion dataframe for week 1
for(i in 1:37){
  assign(paste0("bcr", i, "_mig_prop_week_1"), subset(migratory_proportion_df_week_1, row.names(migratory_proportion_df_week_1) %in% 
                                                        get(paste0("bcr", i, "rows"))))
}

#export dataframes containing migratory proportion for each community, subsetted by BCR
for(i in 1:37){
  if(nrow(get(paste0("bcr", i, "_mig_prop_week_1"))) > 0){ #passes BCRs not represented in week 1 communities
    write.csv(get(paste0("bcr", i, "_mig_prop_week_1")), paste0("output/bcr_community_data/1/bcr", i, "_mig_prop_week_1.csv"))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}

##### VISUALIZE SPATIAL DISTRIBUTION OF MIGRATORY PROPORTION #####

#extract row names defining communities
migratory_proportion_df_week_1.rows <- as.numeric(migratory_proportion_df_week_1$community)

#add coordinates to migratory proportion data frame
migratory_proportion_df_week_1$x <- rep(NA, nrow(migratory_proportion_df_week_1))
migratory_proportion_df_week_1$y <- rep(NA, nrow(migratory_proportion_df_week_1))

#populate coordinates (warning: can take a long time)
for(i in 1:nrow(migratory_proportion_df_week_1)){
  migratory_proportion_df_week_1$x[i] <- xyFromCell(wewpew_1_agg, migratory_proportion_df_week_1.rows[i])[,1]
  migratory_proportion_df_week_1$y[i] <- xyFromCell(wewpew_1_agg, migratory_proportion_df_week_1.rows[i])[,2]
}

#generate gg figure to plot spatial distributions
countries <- ne_countries(returnclass = "sf") %>% 
  st_geometry() %>% 
  st_transform(crs)
country_map <- ggplot()+
  geom_sf(data = countries) +
  coord_sf(xlim=c(-2750000, 3000000), ylim=c(-2000000, 4250000))

#generate map of spatial distribution of migratory proportion in week 1
pdf(file = 'output/figures/community_migratory_proportion_week_1.pdf')
country_map + 
  geom_point(data = migratory_proportion_df_week_1, alpha = 0.1, size = 0.1, aes(x = x, y = y, color = proportion)) +
  scale_color_gradient2("Mig. Prop.", low = "red", high = "blue4",
                        mid = "gray90",
                        midpoint = 0.5) +
  ggtitle("Jan 4") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "")
dev.off()

##### GENERATE COMBINED NRI AND NTI FILES ENCOMPASSING ALL WEEKS #####

#import each NRI and NTI table
for(i in 1:52){
  assign(paste0("week_", i, "_NRI"), read.csv(paste0("output/tables/", i, "/passerine.ses.mpd.", i, ".csv"), row.names=1))
  assign(paste0("week_", i, "_NTI"), read.csv(paste0("output/tables/", i, "/passerine.ses.mntd.", i, ".csv"), row.names=1))
}

## process NRI data
#create a list of all NRI dataframes 
NRIlist <- list(week_1_NRI, week_2_NRI, week_3_NRI, week_4_NRI, week_5_NRI, week_6_NRI, week_7_NRI, week_8_NRI,
                week_9_NRI, week_10_NRI, week_11_NRI, week_12_NRI, week_13_NRI, week_14_NRI, week_15_NRI, week_16_NRI,
                week_17_NRI, week_18_NRI, week_19_NRI, week_20_NRI, week_21_NRI, week_22_NRI, week_23_NRI, week_24_NRI,
                week_25_NRI, week_26_NRI, week_27_NRI, week_28_NRI, week_29_NRI, week_30_NRI, week_31_NRI, week_32_NRI,
                week_33_NRI, week_34_NRI, week_35_NRI, week_36_NRI, week_37_NRI, week_38_NRI, week_39_NRI, week_40_NRI,
                week_41_NRI, week_42_NRI, week_43_NRI, week_44_NRI, week_45_NRI, week_46_NRI, week_47_NRI, week_48_NRI,
                week_49_NRI, week_50_NRI, week_51_NRI, week_52_NRI)

#create a list of all rows for each dataframe
for(i in 1:length(NRIlist)){
  assign(paste0("week_", i, "_NRI_rows"), as.numeric(row.names(get(paste0("week_", i, "_NRI")))))
}
NRIlistrows <- list(week_1_NRI_rows, week_2_NRI_rows, week_3_NRI_rows, week_4_NRI_rows, week_5_NRI_rows, week_6_NRI_rows, week_7_NRI_rows, week_8_NRI_rows,
                    week_9_NRI_rows, week_10_NRI_rows, week_11_NRI_rows, week_12_NRI_rows, week_13_NRI_rows, week_14_NRI_rows, week_15_NRI_rows, week_16_NRI_rows,
                    week_17_NRI_rows, week_18_NRI_rows, week_19_NRI_rows, week_20_NRI_rows, week_21_NRI_rows, week_22_NRI_rows, week_23_NRI_rows, week_24_NRI_rows,
                    week_25_NRI_rows, week_26_NRI_rows, week_27_NRI_rows, week_28_NRI_rows, week_29_NRI_rows, week_30_NRI_rows, week_31_NRI_rows, week_32_NRI_rows,
                    week_33_NRI_rows, week_34_NRI_rows, week_35_NRI_rows, week_36_NRI_rows, week_37_NRI_rows, week_38_NRI_rows, week_39_NRI_rows, week_40_NRI_rows,
                    week_41_NRI_rows, week_42_NRI_rows, week_43_NRI_rows, week_44_NRI_rows, week_45_NRI_rows, week_46_NRI_rows, week_47_NRI_rows, week_48_NRI_rows,
                    week_49_NRI_rows, week_50_NRI_rows, week_51_NRI_rows, week_52_NRI_rows)

#add empty rows that will need to be added
for(i in 1:length(NRIlist)){
  NRIlist[[i]]$x <- rep(NA, nrow(NRIlist[[i]]))
  NRIlist[[i]]$y <- rep(NA, nrow(NRIlist[[i]]))
  NRIlist[[i]]$NRI <- rep(NA, nrow(NRIlist[[i]]))
  NRIlist[[i]]$BCR <- rep(NA, nrow(NRIlist[[i]]))
}

#load species presence file encompassing all communities
wewpew_26_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 26), us), crs, method = "near"), fact=2)

#add coords and NRI to each subsampled community
for(i in 1:length(NRIlist)){
  for(j in 1:nrow(NRIlist[[i]])){
    NRIlist[[i]]$x[j] <- xyFromCell(wewpew_26_agg, NRIlistrows[[i]][j])[,1]
    NRIlist[[i]]$y[j] <- xyFromCell(wewpew_26_agg, NRIlistrows[[i]][j])[,2]
  }
  NRIlist[[i]]$NRI <- -1*NRIlist[[i]]$mpd.obs.z
}

#assign row names to BCRs
for(i in 1:37){
  assign(paste0("bcr",i), allbcrs[allbcrs$BCR == i,] %>%
           st_transform(crs = crs) %>%
           st_geometry())
  assign(paste0("bcr", i, "spatial"), as_Spatial(get(paste0("bcr", i))))
  assign(paste0("wewpew_26_bcr_", i), rast(mask(raster(wewpew_26_agg), get(paste0("bcr", i, "spatial")))))
  assign(paste0("wewpew_26_bcr_", i, "_df"), as.data.frame(get(paste0("wewpew_26_bcr_", i))))
  assign(paste0("bcr", i, "rows"), row.names(get(paste0("wewpew_26_bcr_", i, "_df"))))
}

bcr2rowsdf <- data.frame(rows=bcr2rows, BCR=2)

for(i in 4:37){
  assign(paste0("bcr", i, "rowsdf"), data.frame(rows=get(paste0("bcr", i, "rows")), BCR=i))
}

#combine all BCR rows into one dataframe
allbcrrows <- bind_rows(bcr2rowsdf, bcr4rowsdf, bcr5rowsdf, bcr6rowsdf, bcr7rowsdf, bcr8rowsdf, 
                        bcr9rowsdf, bcr10rowsdf, bcr11rowsdf, bcr12rowsdf, bcr13rowsdf, bcr14rowsdf,
                        bcr15rowsdf, bcr16rowsdf, bcr17rowsdf, bcr18rowsdf, bcr19rowsdf, bcr20rowsdf, 
                        bcr21rowsdf, bcr22rowsdf, bcr23rowsdf, bcr24rowsdf, bcr25rowsdf, bcr26rowsdf,
                        bcr27rowsdf, bcr28rowsdf, bcr29rowsdf, bcr30rowsdf, bcr31rowsdf, bcr32rowsdf,
                        bcr33rowsdf, bcr34rowsdf, bcr35rowsdf, bcr36rowsdf, bcr37rowsdf)

#apply BCR values to NRI list (warning! takes multiple hours)
for(i in 1:length(NRIlist)){
  for(j in 1:nrow(NRIlist[[i]])){
    if((length(allbcrrows[allbcrrows$rows == rownames((NRIlist[[i]])[j,]),"BCR"]) == 1) == TRUE){
      NRIlist[[i]][j,12] <- allbcrrows[allbcrrows$rows == rownames((NRIlist[[i]])[j,]),"BCR"]
    } else {
      print("Pass")
    }
  }
}

#combine into one dataframe
NRIlistdf <- do.call(rbind, unname(NRIlist))

#convert x and y values to lat and long (for downstream analyses with same files)
NRIlistdf$latitude <- sf_project(from = crs, to = st_crs(4326), NRIlistdf[,c("x","y")])[,2]
NRIlistdf$longitude <- sf_project(from = crs, to = st_crs(4326), NRIlistdf[,c("x","y")])[,1]

## process NTI data
#create a list of all NTI dataframes 
NTIlist <- list(week_1_NTI, week_2_NTI, week_3_NTI, week_4_NTI, week_5_NTI, week_6_NTI, week_7_NTI, week_8_NTI,
                week_9_NTI, week_10_NTI, week_11_NTI, week_12_NTI, week_13_NTI, week_14_NTI, week_15_NTI, week_16_NTI,
                week_17_NTI, week_18_NTI, week_19_NTI, week_20_NTI, week_21_NTI, week_22_NTI, week_23_NTI, week_24_NTI,
                week_25_NTI, week_26_NTI, week_27_NTI, week_28_NTI, week_29_NTI, week_30_NTI, week_31_NTI, week_32_NTI,
                week_33_NTI, week_34_NTI, week_35_NTI, week_36_NTI, week_37_NTI, week_38_NTI, week_39_NTI, week_40_NTI,
                week_41_NTI, week_42_NTI, week_43_NTI, week_44_NTI, week_45_NTI, week_46_NTI, week_47_NTI, week_48_NTI,
                week_49_NTI, week_50_NTI, week_51_NTI, week_52_NTI)

#create a list of all rows for each dataframe
for(i in 1:length(NTIlist)){
  assign(paste0("week_", i, "_NTI_rows"), as.numeric(row.names(get(paste0("week_", i, "_NTI")))))
}
NTIlistrows <- list(week_1_NTI_rows, week_2_NTI_rows, week_3_NTI_rows, week_4_NTI_rows, week_5_NTI_rows, week_6_NTI_rows, week_7_NTI_rows, week_8_NTI_rows,
                    week_9_NTI_rows, week_10_NTI_rows, week_11_NTI_rows, week_12_NTI_rows, week_13_NTI_rows, week_14_NTI_rows, week_15_NTI_rows, week_16_NTI_rows,
                    week_17_NTI_rows, week_18_NTI_rows, week_19_NTI_rows, week_20_NTI_rows, week_21_NTI_rows, week_22_NTI_rows, week_23_NTI_rows, week_24_NTI_rows,
                    week_25_NTI_rows, week_26_NTI_rows, week_27_NTI_rows, week_28_NTI_rows, week_29_NTI_rows, week_30_NTI_rows, week_31_NTI_rows, week_32_NTI_rows,
                    week_33_NTI_rows, week_34_NTI_rows, week_35_NTI_rows, week_36_NTI_rows, week_37_NTI_rows, week_38_NTI_rows, week_39_NTI_rows, week_40_NTI_rows,
                    week_41_NTI_rows, week_42_NTI_rows, week_43_NTI_rows, week_44_NTI_rows, week_45_NTI_rows, week_46_NTI_rows, week_47_NTI_rows, week_48_NTI_rows,
                    week_49_NTI_rows, week_50_NTI_rows, week_51_NTI_rows, week_52_NTI_rows)

#add empty rows that will need to be added
for(i in 1:length(NTIlist)){
  NTIlist[[i]]$x <- rep(NA, nrow(NTIlist[[i]]))
  NTIlist[[i]]$y <- rep(NA, nrow(NTIlist[[i]]))
  NTIlist[[i]]$NTI <- rep(NA, nrow(NTIlist[[i]]))
  NTIlist[[i]]$BCR <- rep(NA, nrow(NTIlist[[i]]))
}

#add coords and NTI to each subsampled community
for(i in 1:length(NTIlist)){
  for(j in 1:nrow(NTIlist[[i]])){
    NTIlist[[i]]$x[j] <- xyFromCell(wewpew_26_agg, NTIlistrows[[i]][j])[,1]
    NTIlist[[i]]$y[j] <- xyFromCell(wewpew_26_agg, NTIlistrows[[i]][j])[,2]
  }
  NTIlist[[i]]$NTI <- -1*NTIlist[[i]]$mntd.obs.z
}

#apply BCR values to NTI list (warning! takes takes multiple hours)
for(i in 1:length(NTIlist)){
  for(j in 1:nrow(NTIlist[[i]])){
    if((length(allbcrrows[allbcrrows$rows == rownames((NTIlist[[i]])[j,]),"BCR"]) == 1) == TRUE){
      NTIlist[[i]][j,12] <- allbcrrows[allbcrrows$rows == rownames((NTIlist[[i]])[j,]),"BCR"]
    } else {
      print("Pass")
    }
  }
}

#combine into one dataframe
NTIlistdf <- do.call(rbind, unname(NTIlist))

#convert x and y values to lat and long (for downstream analyses with same files)
NTIlistdf$latitude <- sf_project(from = crs, to = st_crs(4326), NTIlistdf[,c("x","y")])[,2]
NTIlistdf$longitude <- sf_project(from = crs, to = st_crs(4326), NTIlistdf[,c("x","y")])[,1]

##### ADD MIGRATORY PROPORTION TO ALL WEEKS (CONTINETNAL POOLS) #####

#add field for migratory proportion
NRIlistdf$proportion <-  rep(NA, nrow(NRIlistdf))
NTIlistdf$proportion <-  rep(NA, nrow(NTIlistdf))

for(i in 1:52){ #loading in the files takes about a half hour
  assign(paste0("migratory_proportion_df_week", i), read.csv(paste0("output/community_data/migratory_proportion_week_", i, ".csv")))
} 

#combine migratory proportion files into list
all_mig_props <- list(migratory_proportion_df_week1, migratory_proportion_df_week2, migratory_proportion_df_week3, migratory_proportion_df_week4,
                      migratory_proportion_df_week5, migratory_proportion_df_week6, migratory_proportion_df_week7, migratory_proportion_df_week8,
                      migratory_proportion_df_week9, migratory_proportion_df_week10, migratory_proportion_df_week11, migratory_proportion_df_week12,
                      migratory_proportion_df_week13, migratory_proportion_df_week14, migratory_proportion_df_week15, migratory_proportion_df_week16,
                      migratory_proportion_df_week17, migratory_proportion_df_week18, migratory_proportion_df_week19, migratory_proportion_df_week20,
                      migratory_proportion_df_week21, migratory_proportion_df_week22, migratory_proportion_df_week23, migratory_proportion_df_week24,
                      migratory_proportion_df_week25, migratory_proportion_df_week26, migratory_proportion_df_week27, migratory_proportion_df_week28,
                      migratory_proportion_df_week29, migratory_proportion_df_week30, migratory_proportion_df_week31, migratory_proportion_df_week32,
                      migratory_proportion_df_week33, migratory_proportion_df_week34, migratory_proportion_df_week35, migratory_proportion_df_week36,
                      migratory_proportion_df_week37, migratory_proportion_df_week38, migratory_proportion_df_week39, migratory_proportion_df_week40,
                      migratory_proportion_df_week41, migratory_proportion_df_week42, migratory_proportion_df_week43, migratory_proportion_df_week44,
                      migratory_proportion_df_week45, migratory_proportion_df_week46, migratory_proportion_df_week47, migratory_proportion_df_week48,
                      migratory_proportion_df_week49, migratory_proportion_df_week50, migratory_proportion_df_week51, migratory_proportion_df_week52)

#add week column to dataframes
#NRI
NRIlistdf$week <- rep(NA, nrow(NRIlistdf))
NRIlistdf[c(1:25000),"week"] <- 1
for(i in 1:51){
  NRIlistdf[(25000*i)+c(1:25000),"week"] <- i+1
}

#NTI
NTIlistdf$week <- rep(NA, nrow(NTIlistdf))
NTIlistdf[c(1:25000),"week"] <- 1
for(i in 1:51){
  NTIlistdf[(25000*i)+c(1:25000),"week"] <- i+1
}

#populate migratory proportion into NRI data through complete annual cycle
#warning! extremely intensive computationally (~3 days)
for(i in 1:nrow(NRIlistdf)){ 
  NRIlistdf[i,"proportion"] <- all_mig_props[[NRIlistdf[i,"week"]]][all_mig_props[[NRIlistdf[i,"week"]]]$community == NRIlistdf[i,"X"],"proportion"]
}
#populate migratory proportion into NTI data through complete annual cycle
#warning! extremely intensive computationally (~3 days)
for(i in 1:nrow(NTIlistdf)){ 
  NTIlistdf[i,"proportion"] <- all_mig_props[[NTIlistdf[i,"week"]]][all_mig_props[[NTIlistdf[i,"week"]]]$community == NTIlistdf[i,"X"],"proportion"]
}


#export dataframes
write.csv(NRIlistdf, "output/tables/all_NRI_data.csv")
write.csv(NTIlistdf, "output/tables/all_NTI_data.csv")

#reload dataframes if already processed
NRIlistdf <- read.csv("output/tables/all_NRI_data.csv")
NTIlistdf <- read.csv("output/tables/all_NTI_data.csv")

##### ADD MIGRATORY PROPORTION TO ALL WEEKS (ECOREGIONAL POOLS) #####

#import each combined NTI and NRI table
for(i in 1:52){
  assign(paste0("combined_week_", i, "_NRI"), read.csv(paste0("output/bcr_tables/", i, "/combined_NRI_", i, ".csv")))
  assign(paste0("combined_week_", i, "_NTI"), read.csv(paste0("output/bcr_tables/", i, "/combined_NTI_", i, ".csv")))
}

#create list of NRI data for all weeks
all_combined_week_NRIs <- list(combined_week_1_NRI, combined_week_2_NRI, combined_week_3_NRI, combined_week_4_NRI,
                               combined_week_5_NRI, combined_week_6_NRI, combined_week_7_NRI, combined_week_8_NRI,
                               combined_week_9_NRI, combined_week_10_NRI, combined_week_11_NRI, combined_week_12_NRI,
                               combined_week_13_NRI, combined_week_14_NRI, combined_week_15_NRI, combined_week_16_NRI,
                               combined_week_17_NRI, combined_week_18_NRI, combined_week_19_NRI, combined_week_20_NRI,
                               combined_week_21_NRI, combined_week_22_NRI, combined_week_23_NRI, combined_week_24_NRI,
                               combined_week_25_NRI, combined_week_26_NRI, combined_week_27_NRI, combined_week_28_NRI,
                               combined_week_29_NRI, combined_week_30_NRI, combined_week_31_NRI, combined_week_32_NRI,
                               combined_week_33_NRI, combined_week_34_NRI, combined_week_35_NRI, combined_week_36_NRI,
                               combined_week_37_NRI, combined_week_38_NRI, combined_week_39_NRI, combined_week_40_NRI,
                               combined_week_41_NRI, combined_week_42_NRI, combined_week_43_NRI, combined_week_44_NRI,
                               combined_week_45_NRI, combined_week_46_NRI, combined_week_47_NRI, combined_week_48_NRI,
                               combined_week_49_NRI, combined_week_50_NRI, combined_week_51_NRI, combined_week_52_NRI)

#create list of NTI data for all weeks
all_combined_week_NTIs <- list(combined_week_1_NTI, combined_week_2_NTI, combined_week_3_NTI, combined_week_4_NTI,
                               combined_week_5_NTI, combined_week_6_NTI, combined_week_7_NTI, combined_week_8_NTI,
                               combined_week_9_NTI, combined_week_10_NTI, combined_week_11_NTI, combined_week_12_NTI,
                               combined_week_13_NTI, combined_week_14_NTI, combined_week_15_NTI, combined_week_16_NTI,
                               combined_week_17_NTI, combined_week_18_NTI, combined_week_19_NTI, combined_week_20_NTI,
                               combined_week_21_NTI, combined_week_22_NTI, combined_week_23_NTI, combined_week_24_NTI,
                               combined_week_25_NTI, combined_week_26_NTI, combined_week_27_NTI, combined_week_28_NTI,
                               combined_week_29_NTI, combined_week_30_NTI, combined_week_31_NTI, combined_week_32_NTI,
                               combined_week_33_NTI, combined_week_34_NTI, combined_week_35_NTI, combined_week_36_NTI,
                               combined_week_37_NTI, combined_week_38_NTI, combined_week_39_NTI, combined_week_40_NTI,
                               combined_week_41_NTI, combined_week_42_NTI, combined_week_43_NTI, combined_week_44_NTI,
                               combined_week_45_NTI, combined_week_46_NTI, combined_week_47_NTI, combined_week_48_NTI,
                               combined_week_49_NTI, combined_week_50_NTI, combined_week_51_NTI, combined_week_52_NTI)

#turn into data frames
all_combined_week_NRIs_df <- data.table::rbindlist(all_combined_week_NRIs, use.names=TRUE, idcol=TRUE, fill=TRUE) %>% 
  as.data.frame()
all_combined_week_NTIs_df <- data.table::rbindlist(all_combined_week_NTIs, use.names=TRUE, idcol=TRUE, fill=TRUE) %>% 
  as.data.frame()

#convert x and y values to lat and long (for downstream analyses with same files)
all_combined_week_NRIs_df$latitude <- sf_project(from = crs, to = st_crs(4326), all_combined_week_NRIs_df[,c("x","y")])[,2]
all_combined_week_NRIs_df$longitude <- sf_project(from = crs, to = st_crs(4326), all_combined_week_NRIs_df[,c("x","y")])[,1]
all_combined_week_NTIs_df$latitude <- sf_project(from = crs, to = st_crs(4326), all_combined_week_NTIs_df[,c("x","y")])[,2]
all_combined_week_NTIs_df$longitude <- sf_project(from = crs, to = st_crs(4326), all_combined_week_NTIs_df[,c("x","y")])[,1]

## process for NRI data
#add field for migratory proportion
all_combined_week_NRIs_df$proportion <- rep(NA, nrow(all_combined_week_NRIs_df))

#populate migratory proportion into NRI data through complete annual cycle
#warning! extremely intensive computationally (~3 days)
for(i in 1:nrow(all_combined_week_NRIs_df)){ 
  all_combined_week_NRIs_df[i,"proportion"] <- all_mig_props[[all_combined_week_NRIs_df[i,".id"]]][all_mig_props[[all_combined_week_NRIs_df[i,".id"]]]$community == all_combined_week_NRIs_df[i,"X"],"proportion"]
}

#export updated NRI data file containing migratory proportion
write.csv(all_combined_week_NRIs_df, "output/tables/combined_NRI_df.csv")

#reload dataframe if already processed
all_combined_week_NRIs_df <- read.csv("output/tables/combined_NRI_df.csv")

## replicate process with NTI data
#add field for migratory proportion
all_combined_week_NTIs_df$proportion <- rep(NA, nrow(all_combined_week_NTIs_df))

#populate migratory proportion into NRI data through complete annual cycle
#warning! extremely intensive computationally (~3 days)
for(i in 1:nrow(all_combined_week_NTIs_df)){ 
  all_combined_week_NTIs_df[i,"proportion"] <- all_mig_props[[all_combined_week_NTIs_df[i,".id"]]][all_mig_props[[all_combined_week_NTIs_df[i,".id"]]]$community == all_combined_week_NTIs_df[i,"X"],"proportion"]
}

#export updated NTI data file containing migratory proportion
write.csv(all_combined_week_NTIs_df, "output/tables/combined_NTI_df.csv")

#reload dataframe if already processed
all_combined_week_NTIs_df <- read.csv("output/tables/combined_NTI_df.csv")
