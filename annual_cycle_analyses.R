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
library(rgdal)
library(raster)
library(devtools)
library(data.table)
#devtools::install_github("ropensci/rnaturalearthhires")
extract <- terra::extract
project <- terra::project #possible source of confusion regarding rgdal for later

# Check for needed packages, and install the missing ones
required.libraries <- c("ape", "picante", 
                        "pez", "phytools",
                        "vegan", "adephylo", 
                        "phylobase", "geiger", 
                        "mvMORPH", "OUwie", 
                        "hisse", "BAMMtools",
                        "phylosignal", "Biostrings",
                        "devtools","ggplot2", 
                        "kableExtra", "betapart", "gridExtra",
                        "reshape2")

# Pure intro stuff (load to help set parameters and define common ebirdst terms)

# download a simplified example dataset for Yellow-bellied Sapsucker in Michigan
path <- ebirdst_download(species = "example_data")

# for non-example data use the species code or name instead of "example_data"
path <- get_species_path("example_data")

occ_lr <- load_raster(path, product = "occurrence", resolution = "lr")
occ_lr

parse_raster_dates(occ_lr)

# weekly, low res, abundance confidence intervals
abd_lower <- load_raster(path, product = "abundance", metric = "lower", 
                         resolution = "lr")
abd_upper <- load_raster(path, product = "abundance", metric = "upper", 
                         resolution = "lr")


# seasonal, low res, mean relative abundance
abd_seasonal_mean <- load_raster(path, product = "abundance", 
                                 period = "seasonal", metric = "mean", 
                                 resolution = "lr")

# season that each layer corresponds to
names(abd_seasonal_mean)

# just the breeding season layer
abd_seasonal_mean[["breeding"]]

# seasonal, low res, max occurrence
occ_seasonal_max <- load_raster(path, product = "occurrence", 
                                period = "seasonal", metric = "max", 
                                resolution = "lr")

# full year, low res, maximum relative abundance
abd_fy_max <- load_raster(path, product = "abundance", 
                          period = "full-year", metric = "max", 
                          resolution = "lr")


# Working with rasters module

# load seasonal mean relative abundance at low res
abd_seasonal <- load_raster(path, 
                            product = "abundance", 
                            period = "seasonal",
                            metric = "mean",
                            resolution = "lr")

abd_breeding <- abd_seasonal[["breeding"]]

# boundary polygon for michigan
mi <- ne_states(iso_a2 = "US", returnclass = "sf") %>% 
  filter(postal == "MI") %>% 
  # project to same coordinate reference system as the raster data
  st_transform(st_crs(abd_seasonal))

# crop data to michigan
abd_breeding_mi <- crop(abd_breeding, mi)

# load the mapping parameters
fac_parameters <- load_fac_map_parameters(path)
crs <- fac_parameters$custom_projection

# transform to the custom projection using nearest neighbor resampling
abd_projected <- project(abd_breeding_mi, crs, method = "near")

# quantiles of non-zero values
v <- values(abd_projected)
v <- v[!is.na(v) & v > 0]
bins <- quantile(v, seq(0, 1, by = 0.1))
# add a bin for 0
bins <- c(0, bins)

# status and trends palette
pal <- abundance_palette(length(bins) - 2)
# add a color for zero
pal <- c("#e6e6e6", pal)

# natural earth boundaries
countries <- ne_countries(returnclass = "sf") %>% 
  st_geometry() %>% 
  st_transform(crs)
states <- ne_states(iso_a2 = "US", returnclass = "sf") %>% 
  st_geometry() %>% 
  st_transform(crs)


# define the map extent with the michigan polygon
mi_ext <- mi %>% 
  st_geometry() %>% 
  st_transform(crs)

labels <- quantile(bins, c(0, 0.5, 1))
label_breaks <- seq(0, 1, length.out = length(bins))

abd_median <- load_raster(path, product = "abundance", 
                          metric = "median", resolution = "lr")
abd_lower <- load_raster(path, product = "abundance",
                         metric = "lower", resolution = "lr")
abd_upper <- load_raster(path, product = "abundance", 
                         metric = "upper", resolution = "lr")

#set some geographic parameters

ovenbi1_abd_seasonal_mean <- load_raster("data/ebirdst/2021/ovenbi1", product = "abundance", 
                                         period = "seasonal", metric = "mean", 
                                         resolution = "hr")


us <- ne_countries(country = "united states of america", type = "countries", returnclass = "sf") %>% 
  st_transform(st_crs(ovenbi1_abd_seasonal_mean))

#import each combined NTI and NRI table
for(i in 1:52){
  assign(paste0("combined_week_", i, "_NRI"), read.csv(paste0("output/bcr_tables/", i, "/combined_NRI_", i, ".csv")))
  assign(paste0("combined_week_", i, "_NTI"), read.csv(paste0("output/bcr_tables/", i, "/combined_NTI_", i, ".csv")))
}



#create list of all combined week NRIs
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



#replicate with list of all combined week NTIs
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

#create table of mean NRI values
BCRs <- c(5:37)
weeks <- c(1:52)
mean_NRI_values_matrix <- matrix(NA, ncol = length(weeks), nrow = length(BCRs))
col_names <- as.character(weeks)

# Convert the matrix to a data frame with column names
mean_NRI_values <- data.frame(matrix = mean_NRI_values_matrix)
colnames(mean_NRI_values) <- col_names
rownames(mean_NRI_values) <- BCRs

#populate the mean NRI values table
for(i in BCRs){ 
  for(j in weeks){
    mean_NRI_values[i-4, j] <- mean(na.omit(all_combined_week_NRIs[[j]][all_combined_week_NRIs[[j]]$BCR == i, "NRI"]))
  }
}

mean_NRI_values$BCR <- BCRs


#create table of mean NTI values
mean_NTI_values_matrix <- matrix(NA, ncol = length(weeks), nrow = length(BCRs))

# Convert the matrix to a data frame with column names
mean_NTI_values <- data.frame(matrix = mean_NTI_values_matrix)
colnames(mean_NTI_values) <- col_names
rownames(mean_NTI_values) <- BCRs

#populate the mean NTI values table
for(i in BCRs){ 
  for(j in weeks){
    mean_NTI_values[i-4, j] <- mean(na.omit(all_combined_week_NTIs[[j]][all_combined_week_NTIs[[j]]$BCR == i, "NTI"]))
  }
}

mean_NTI_values$BCR <- BCRs


#switch off spherical geometry
sf_use_s2(FALSE)

# get BCR latitude values
#create BCR centroids files
allbcrs <- read_sf("data/BCR_Terrestrial/BCR_Terrestrial_master.shp") #read in BCR files
for(i in BCRs){
  assign(paste0("bcr_centroids_", i), st_centroid(allbcrs[allbcrs$BCR == i,"geometry"]) %>% 
           st_as_sf(wkt = "a") %>% 
           st_coordinates())
}

#get the mean of the centroids
mean(bcr_centroids_5[,2]) #50.6826
mean(bcr_centroids_6[,2]) #58.50326
mean(bcr_centroids_7[,2]) #57.83681
mean(bcr_centroids_8[,2]) #54.43946
mean(bcr_centroids_9[,2]) #43.24103
mean(bcr_centroids_10[,2]) #46.19079
mean(bcr_centroids_11[,2]) #47.40598
mean(bcr_centroids_12[,2]) #47.95409
mean(bcr_centroids_13[,2]) #43.65609
mean(bcr_centroids_14[,2]) #44.72067
mean(bcr_centroids_15[,2]) #38.88149
mean(bcr_centroids_16[,2]) #38.52856
mean(bcr_centroids_17[,2]) #44.91282
mean(bcr_centroids_18[,2]) #38.69052
mean(bcr_centroids_19[,2]) #38.43332
mean(bcr_centroids_20[,2]) #30.43795
mean(bcr_centroids_21[,2]) #34.65978
mean(bcr_centroids_22[,2]) #40.8388
mean(bcr_centroids_23[,2]) #42.93611
mean(bcr_centroids_24[,2]) #37.05901
mean(bcr_centroids_25[,2]) #33.09216
mean(bcr_centroids_26[,2]) #35.14722
mean(bcr_centroids_27[,2]) #33.5757
mean(bcr_centroids_28[,2]) #37.8537
mean(bcr_centroids_29[,2]) #37.49024
mean(bcr_centroids_30[,2]) #40.89619
mean(bcr_centroids_31[,2]) #27.68041
mean(bcr_centroids_32[,2]) #34.09171
mean(bcr_centroids_33[,2]) #32.76825
mean(bcr_centroids_34[,2]) #25.70883
mean(bcr_centroids_35[,2]) #26.28559
mean(bcr_centroids_36[,2]) #26.81768
mean(bcr_centroids_37[,2]) #28.57299

#add latitudes to dataframes
mean_latitudes <- c(50.6826, 58.50326, 57.83681, 54.43946, 43.24103, 46.19079,
                    47.40598, 47.95409, 43.65609, 44.72067, 38.88149, 38.52856,
                    44.91282, 38.69052, 38.43332, 30.43795, 34.65978, 40.8388,
                    42.93611, 37.05901, 33.09216, 35.14722, 33.5757, 37.8537,
                    37.49024, 40.89619, 27.68041, 34.09171, 32.76825, 25.70883,
                    26.28559, 26.81768, 28.57299)

mean_NRI_values$latitude <- mean_latitudes
mean_NTI_values$latitude <- mean_latitudes



# attempt visualizations

#NRI
mean_NRI_values_visuals <- reshape2::melt(mean_NRI_values, id.vars=c("BCR", "latitude"))
mean_NRI_values_visuals

ggplot(mean_NRI_values_visuals, aes(x = variable, y = value)) + geom_smooth(aes(group = BCR, colour=latitude), se=F) + 
  scale_colour_viridis(direction = 1) + 
  ggtitle("Annual NRI Distribution by BCR") + 
  labs(x = "Week", y = "NRI")



#NTI
mean_NTI_values_visuals <- reshape2::melt(mean_NTI_values, id.vars=c("BCR", "latitude"))
mean_NTI_values_visuals

#visualize
ggplot(mean_NTI_values_visuals, aes(x = variable, y = value)) + geom_smooth(aes(group = BCR, colour=latitude), se=F) + 
  scale_colour_viridis(direction = 1) + 
  ggtitle("Annual NTI Distribution by BCR") + 
  labs(x = "Week", y = "NTI")



### Get weekly BCR species diversity
#split it up seasonally so that R can handle the amount of data without crashing

#weeks 1-10
#for(j in 1:10){
#  for(i in 1:37){
#    assign(paste0("community_", j, "_bcr_", i), data.table::fread(file = paste0("output/bcr_community_data/", j, "/community_", j, "_bcr_", i, ".csv")))
#  }
#}

#weeks_1_to_10_diversity <- matrix(NA, ncol = 10, nrow = 37)
#rownames(weeks_1_to_10_diversity) <- c(1:37)
#colnames(weeks_1_to_10_diversity) <- c(1:10)


#for(j in 1:ncol(weeks_1_to_10_diversity)){
#  for(i in 1:nrow(weeks_1_to_10_diversity)){
#    weeks_1_to_10_diversity[i,j] <- ncol(get(paste0("community_", j, "_bcr_", i)))
#  }
#}

#write.csv(weeks_1_to_10_diversity, "output/tables/weeks_1_to_10_bcr_diversity.csv")
weeks_1_to_10_diversity <- read.csv("output/tables/weeks_1_to_10_bcr_diversity.csv", row.names=1)


#weeks 11-20
#for(j in 11:20){
#  for(i in 1:37){
#    assign(paste0("community_", j, "_bcr_", i), data.table::fread(file = paste0("output/bcr_community_data/", j, "/community_", j, "_bcr_", i, ".csv")))
#  }
#}

#weeks_11_to_20_diversity <- matrix(NA, ncol = 10, nrow = 37)
#rownames(weeks_11_to_20_diversity) <- c(1:37)
#colnames(weeks_11_to_20_diversity) <- c(11:20)


#for(j in 11:20){
#  for(i in 1:nrow(weeks_11_to_20_diversity)){
#    weeks_11_to_20_diversity[i,j-10] <- ncol(get(paste0("community_", j, "_bcr_", i)))
#  }
#}

#write.csv(weeks_11_to_20_diversity, "output/tables/weeks_11_to_20_bcr_diversity.csv")
weeks_11_to_20_diversity <- read.csv("output/tables/weeks_11_to_20_bcr_diversity.csv", row.names=1)



#weeks 21-30
#for(j in 21:30){
#  for(i in 1:37){
#    assign(paste0("community_", j, "_bcr_", i), data.table::fread(file = paste0("output/bcr_community_data/", j, "/community_", j, "_bcr_", i, ".csv")))
#  }
#}

#weeks_21_to_30_diversity <- matrix(NA, ncol = 10, nrow = 37)
#rownames(weeks_21_to_30_diversity) <- c(1:37)
#colnames(weeks_21_to_30_diversity) <- c(21:30)


#for(j in 21:30){
#  for(i in 1:nrow(weeks_21_to_30_diversity)){
#    weeks_21_to_30_diversity[i,j-20] <- ncol(get(paste0("community_", j, "_bcr_", i)))
#  }
#}

#write.csv(weeks_21_to_30_diversity, "output/tables/weeks_21_to_30_bcr_diversity.csv")
weeks_21_to_30_diversity <- read.csv("output/tables/weeks_21_to_30_bcr_diversity.csv", row.names=1)



#weeks 31-40
#for(j in 31:40){
#  for(i in 1:37){
#    assign(paste0("community_", j, "_bcr_", i), data.table::fread(file = paste0("output/bcr_community_data/", j, "/community_", j, "_bcr_", i, ".csv")))
#  }
#}

#weeks_31_to_40_diversity <- matrix(NA, ncol = 10, nrow = 37)
#rownames(weeks_31_to_40_diversity) <- c(1:37)
#colnames(weeks_31_to_40_diversity) <- c(31:40)


#for(j in 31:40){
#  for(i in 1:nrow(weeks_31_to_40_diversity)){
#    weeks_31_to_40_diversity[i,j-30] <- ncol(get(paste0("community_", j, "_bcr_", i)))
#  }
#}

#write.csv(weeks_31_to_40_diversity, "output/tables/weeks_31_to_40_bcr_diversity.csv")
weeks_31_to_40_diversity <- read.csv("output/tables/weeks_31_to_40_bcr_diversity.csv", row.names=1)



#weeks 41-52
#for(j in 41:52){
#  for(i in 1:37){
#    assign(paste0("community_", j, "_bcr_", i), data.table::fread(file = paste0("output/bcr_community_data/", j, "/community_", j, "_bcr_", i, ".csv")))
#  }
#}

#weeks_41_to_52_diversity <- matrix(NA, ncol = 12, nrow = 37)
#rownames(weeks_41_to_52_diversity) <- c(1:37)
#colnames(weeks_41_to_52_diversity) <- c(41:52)


#for(j in 41:52){
#  for(i in 1:nrow(weeks_41_to_52_diversity)){
#    weeks_41_to_52_diversity[i,j-40] <- ncol(get(paste0("community_", j, "_bcr_", i)))
#  }
#}

#write.csv(weeks_41_to_52_diversity, "output/tables/weeks_41_to_52_bcr_diversity.csv")
weeks_41_to_52_diversity <- read.csv("output/tables/weeks_41_to_52_bcr_diversity.csv", row.names=1)
weeks_41_to_52_diversity

#weekly_bcr_diversity <- bind_cols(weeks_1_to_10_diversity, weeks_11_to_20_diversity, weeks_21_to_30_diversity, weeks_31_to_40_diversity, weeks_41_to_52_diversity)
#colnames(weekly_bcr_diversity) <- c(1:52)
#weekly_bcr_diversity[1,26] <- 0
#weekly_bcr_diversity[3,26] <- 0
#write.csv(weekly_bcr_diversity, "output/tables/weekly_bcr_diversity.csv")
weekly_bcr_diversity <- read.csv("output/tables/weekly_bcr_diversity.csv")
weekly_bcr_diversity


weekly_bcr_diversity_visuals <- weekly_bcr_diversity[5:37,2:53]
weekly_bcr_diversity_values <- as.vector(unlist(weekly_bcr_diversity_visuals))
weekly_bcr_diversity_visuals

mean_NRI_values_visuals$richness <- weekly_bcr_diversity_values
mean_NTI_values_visuals$richness <- weekly_bcr_diversity_values


### Get weekly BCR migratory proportion total
#split it up seasonally so that R can handle the amount of data without crashing

#weeks 1-10
# for(j in 1:10){
#   for(i in 1:37){
#     f <- paste0("output/bcr_community_data/", j, "/bcr", i, "_mig_prop_week_", j, ".csv") #check if file exists
#     if(file.exists(f)){
#       assign(paste0("community_", j, "_bcr_", i, "_migratory_proportion"), data.table::fread(f)) #read in files that exist
#     }
#     else {
#       message(paste0("Skipped BCR ", i)) #show skipped BCRs
#     }
#   }
# }
# 
# weeks_1_to_10_mig_prop <- matrix(NA, ncol = 10, nrow = 37)
# rownames(weeks_1_to_10_mig_prop) <- c(1:37)
# colnames(weeks_1_to_10_mig_prop) <- c(1:10)
# 
# for(j in 1:ncol(weeks_1_to_10_mig_prop)){
#   for(i in 1:nrow(weeks_1_to_10_mig_prop)){
#     f <- paste0("output/bcr_community_data/", j, "/bcr", i, "_mig_prop_week_", j, ".csv") #check if file exists
#     if(file.exists(f)) {
#       m <- get(paste0("community_", j, "_bcr_", i, "_migratory_proportion"))
#       m <- mean(m$proportion[!is.infinite(m$proportion)])
#       weeks_1_to_10_mig_prop[i,j] <- m
#     }
#     else {
#       weeks_1_to_10_mig_prop[i,j] <- 0
#     }
#   }
# }

#write.csv(weeks_1_to_10_mig_prop, "output/tables/weeks_1_to_10_mig_prop.csv")
weeks_1_to_10_mig_prop <- read.csv("output/tables/weeks_1_to_10_mig_prop.csv", row.names=1)


#weeks 11-20
# for(j in 11:20){
#   for(i in 1:37){
#     f <- paste0("output/bcr_community_data/", j, "/bcr", i, "_mig_prop_week_", j, ".csv") #check if file exists
#     if(file.exists(f)){
#       assign(paste0("community_", j, "_bcr_", i, "_migratory_proportion"), data.table::fread(f)) #read in files that exist
#     }
#     else {
#       message(paste0("Skipped BCR ", i)) #show skipped BCRs
#     }
#   }
# }
# 
# weeks_11_to_20_mig_prop <- matrix(NA, ncol = 10, nrow = 37)
# rownames(weeks_11_to_20_mig_prop) <- c(1:37)
# colnames(weeks_11_to_20_mig_prop) <- c(11:20)
# 
# for(j in 11:20){
#   for(i in 1:nrow(weeks_11_to_20_mig_prop)){
#     f <- paste0("output/bcr_community_data/", j, "/bcr", i, "_mig_prop_week_", j, ".csv") #check if file exists
#     if(file.exists(f)) {
#       m <- get(paste0("community_", j, "_bcr_", i, "_migratory_proportion"))
#       m <- mean(m$proportion[!is.infinite(m$proportion)])
#       weeks_11_to_20_mig_prop[i,j-10] <- m
#     }
#     else {
#       weeks_11_to_20_mig_prop[i,j-10] <- 0
#     }
#   }
# }
# 
# write.csv(weeks_11_to_20_mig_prop, "output/tables/weeks_11_to_20_mig_prop.csv")
weeks_11_to_20_mig_prop <- read.csv("output/tables/weeks_11_to_20_mig_prop.csv", row.names=1)


#weeks 21-30
# for(j in 21:30){
#   for(i in 1:37){
#     f <- paste0("output/bcr_community_data/", j, "/bcr", i, "_mig_prop_week_", j, ".csv") #check if file exists
#     if(file.exists(f)){
#       assign(paste0("community_", j, "_bcr_", i, "_migratory_proportion"), data.table::fread(f)) #read in files that exist
#     }
#     else {
#       message(paste0("Skipped BCR ", i)) #show skipped BCRs
#     }
#   }
# }
# 
# weeks_21_to_30_mig_prop <- matrix(NA, ncol = 10, nrow = 37)
# rownames(weeks_21_to_30_mig_prop) <- c(1:37)
# colnames(weeks_21_to_30_mig_prop) <- c(21:30)
# 
# for(j in 21:30){
#   for(i in 1:nrow(weeks_21_to_30_mig_prop)){
#     f <- paste0("output/bcr_community_data/", j, "/bcr", i, "_mig_prop_week_", j, ".csv") #check if file exists
#     if(file.exists(f)) {
#       m <- get(paste0("community_", j, "_bcr_", i, "_migratory_proportion"))
#       m <- mean(m$proportion[!is.infinite(m$proportion)])
#       weeks_21_to_30_mig_prop[i,j-20] <- m
#     }
#     else {
#       weeks_21_to_30_mig_prop[i,j-20] <- 0
#     }
#   }
# }
# 
# write.csv(weeks_21_to_30_mig_prop, "output/tables/weeks_21_to_30_mig_prop.csv")
weeks_21_to_30_mig_prop <- read.csv("output/tables/weeks_21_to_30_mig_prop.csv", row.names=1)


#weeks 31-40
# for(j in 31:40){
#   for(i in 1:37){
#     f <- paste0("output/bcr_community_data/", j, "/bcr", i, "_mig_prop_week_", j, ".csv") #check if file exists
#     if(file.exists(f)){
#       assign(paste0("community_", j, "_bcr_", i, "_migratory_proportion"), data.table::fread(f)) #read in files that exist
#     }
#     else {
#       message(paste0("Skipped BCR ", i)) #show skipped BCRs
#     }
#   }
# }
# 
# weeks_31_to_40_mig_prop <- matrix(NA, ncol = 10, nrow = 37)
# rownames(weeks_31_to_40_mig_prop) <- c(1:37)
# colnames(weeks_31_to_40_mig_prop) <- c(31:40)
# 
# for(j in 31:40){
#   for(i in 1:nrow(weeks_31_to_40_mig_prop)){
#     f <- paste0("output/bcr_community_data/", j, "/bcr", i, "_mig_prop_week_", j, ".csv") #check if file exists
#     if(file.exists(f)) {
#       m <- get(paste0("community_", j, "_bcr_", i, "_migratory_proportion"))
#       m <- mean(m$proportion[!is.infinite(m$proportion)])
#       weeks_31_to_40_mig_prop[i,j-30] <- m
#     }
#     else {
#       weeks_31_to_40_mig_prop[i,j-30] <- 0
#     }
#   }
# }
# 
# write.csv(weeks_31_to_40_mig_prop, "output/tables/weeks_31_to_40_mig_prop.csv")
weeks_31_to_40_mig_prop <- read.csv("output/tables/weeks_31_to_40_mig_prop.csv", row.names=1)


#weeks 41-52
# for(j in 41:52){
#   for(i in 1:37){
#     f <- paste0("output/bcr_community_data/", j, "/bcr", i, "_mig_prop_week_", j, ".csv") #check if file exists
#     if(file.exists(f)){
#       assign(paste0("community_", j, "_bcr_", i, "_migratory_proportion"), data.table::fread(f)) #read in files that exist
#     }
#     else {
#       message(paste0("Skipped BCR ", i)) #show skipped BCRs
#     }
#   }
# }
# 
# weeks_41_to_52_mig_prop <- matrix(NA, ncol = 12, nrow = 37)
# rownames(weeks_41_to_52_mig_prop) <- c(1:37)
# colnames(weeks_41_to_52_mig_prop) <- c(41:52)
# 
# for(j in 41:52){
#   for(i in 1:nrow(weeks_41_to_52_mig_prop)){
#     f <- paste0("output/bcr_community_data/", j, "/bcr", i, "_mig_prop_week_", j, ".csv") #check if file exists
#     if(file.exists(f)) {
#       m <- get(paste0("community_", j, "_bcr_", i, "_migratory_proportion"))
#       m <- mean(m$proportion[!is.infinite(m$proportion)])
#       weeks_41_to_52_mig_prop[i,j-40] <- m
#     }
#     else {
#       weeks_41_to_52_mig_prop[i,j-40] <- 0
#     }
#   }
# }
# 
# write.csv(weeks_41_to_52_mig_prop, "output/tables/weeks_41_to_52_mig_prop.csv")
weeks_41_to_52_mig_prop <- read.csv("output/tables/weeks_41_to_52_mig_prop.csv", row.names=1)


# add migratory proportion data and make visualizable

# weekly_bcr_mig_prop <- bind_cols(weeks_1_to_10_mig_prop, weeks_11_to_20_mig_prop, weeks_21_to_30_mig_prop, weeks_31_to_40_mig_prop, weeks_41_to_52_mig_prop)
# colnames(weekly_bcr_mig_prop) <- c(1:52)
# 
# write.csv(weekly_bcr_mig_prop, "output/tables/weekly_bcr_mig_prop.csv")
weekly_bcr_mig_prop <- read.csv("output/tables/weekly_bcr_mig_prop.csv")


weekly_bcr_mig_prop_visuals <- weekly_bcr_mig_prop[5:37,2:53]
weekly_bcr_mig_prop_values <- as.vector(unlist(weekly_bcr_mig_prop_visuals))

mean_NRI_values_visuals$mig_prop <- weekly_bcr_mig_prop_values
mean_NTI_values_visuals$mig_prop <- weekly_bcr_mig_prop_values

hist(mean_NTI_values_visuals$mig_prop)


# attempt to visualize richness relationship

#visualize
ggplot(mean_NTI_values_visuals, aes(x = variable, y = value)) + 
  geom_point(aes(colour = richness)) + 
  geom_smooth(aes(group = BCR, colour = richness), se=F) + 
  scale_colour_viridis(direction = -1) + 
  ggtitle("Annual NTI Distribution by BCR") + 
  labs(x = "Week", y = "NTI")

?geom_smooth()

library(tidytext)
#visualize

### facet wrap of weekly NRI values by BCR, points colored by richness
ggplot(mean_NRI_values_visuals, aes(x = variable, y = value)) + 
  geom_point(aes(colour = richness)) + 
  geom_smooth(aes(group = BCR), se=F) + 
  facet_wrap(~BCR) + 
  scale_colour_viridis(direction = 1) + 
  ggtitle("Annual NRI Distribution by BCR") + 
  labs(x = "Week", y = "NRI")

#overall NRI richness relationship
ggplot(mean_NRI_values_visuals, aes(x = richness, y = value)) + 
  geom_point(aes(colour = BCR)) + 
  scale_colour_viridis(direction = 1) +
  geom_smooth(se=F)

#overall NTI richness relationship
ggplot(mean_NTI_values_visuals, aes(x = richness, y = value)) + 
  geom_point(aes(colour = BCR)) + 
  scale_colour_viridis(direction = 1) +
  geom_smooth(aes(colour = ..y..), se=F)

#overall NRI mig prop relationship
ggplot(mean_NRI_values_visuals, aes(x = mig_prop, y = value)) + 
  geom_point(aes(colour = BCR)) + 
  scale_colour_viridis(direction = 1) +
  geom_smooth(se=F)

#overall NTI mig prop relationship
ggplot(mean_NTI_values_visuals, aes(x = mig_prop, y = value)) + 
  geom_point(aes(colour = latitude)) + 
  scale_colour_viridis(direction = 1) +
  geom_smooth(se=F)


### facet wrap of weekly NTI values by BCR, boreal forest
ggplot(filter(mean_NTI_values_visuals, BCR %in% c(6,7,8,12)), aes(x = variable, y = value)) + 
  geom_point(aes(colour = richness), size=3) + 
  geom_smooth(aes(group = BCR), se=F) + 
  facet_wrap(~BCR) + 
  scale_colour_viridis(direction = 1) + 
  ggtitle("Annual NTI Distribution by BCR -- Boreal Forests") + 
  labs(x = "Week", y = "NTI") +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(limits = c(1:52), breaks = (seq(1,49, by=4)))

### facet wrap of weekly NTI values by BCR, southern US forests
ggplot(filter(mean_NTI_values_visuals, BCR %in% c(25,26,27,31,37)), aes(x = variable, y = value)) + 
  geom_point(aes(colour = richness), size=3) + 
  geom_smooth(aes(group = BCR), se=F) + 
  facet_wrap(~BCR) + 
  scale_colour_viridis(direction = 1) + 
  ggtitle("Annual NTI Distribution by BCR -- Southern Forests") + 
  labs(x = "Week", y = "NTI") +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(limits = c(1:52), breaks = (seq(1,49, by=4)))


### facet wrap of weekly NRI values by BCR, boreal forest
ggplot(filter(mean_NRI_values_visuals, BCR %in% c(6,7,8,12)), aes(x = variable, y = value)) + 
  geom_point(aes(colour = richness), size=3) + 
  geom_smooth(aes(group = BCR), se=F) + 
  facet_wrap(~BCR) + 
  scale_colour_viridis(direction = 1) + 
  ggtitle("Annual NRI Distribution by BCR -- Boreal Forests") + 
  labs(x = "Week", y = "NRI")

### facet wrap of weekly NRI values by BCR, southern US forests
ggplot(filter(mean_NRI_values_visuals, BCR %in% c(25,26,27,31,37)), aes(x = variable, y = value)) + 
  geom_point(aes(colour = richness), size=3) + 
  geom_smooth(aes(group = BCR), se=F) + 
  facet_wrap(~BCR) + 
  scale_colour_viridis(direction = 1) + 
  ggtitle("Annual NRI Distribution by BCR -- Southern Forests") + 
  labs(x = "Week", y = "NRI")

#create NRI, NTI, and richness difference vectors
week1_NRIs <- mean_NRI_values_visuals[mean_NRI_values_visuals$variable == 1,"value"]
week1_NTIs <- mean_NTI_values_visuals[mean_NTI_values_visuals$variable == 1,"value"]
week1_richness <- mean_NRI_values_visuals[mean_NRI_values_visuals$variable == 1,"richness"]
week26_NRIs <- mean_NRI_values_visuals[mean_NRI_values_visuals$variable == 26,"value"]
week26_NTIs <- mean_NTI_values_visuals[mean_NTI_values_visuals$variable == 26,"value"]
week26_richness <- mean_NRI_values_visuals[mean_NRI_values_visuals$variable == 26,"richness"]
NRI_diff <- week26_NRIs - week1_NRIs
NTI_diff <- week26_NTIs - week1_NTIs
richness_diff <- week26_richness - week1_richness

#create mig_prop richness difference vectors
week1_mig_prop <- mean_NRI_values_visuals[mean_NRI_values_visuals$variable == 1,"mig_prop"]
week26_mig_prop <- mean_NRI_values_visuals[mean_NRI_values_visuals$variable == 26,"mig_prop"]
mig_prop_diff <- week26_mig_prop - week1_mig_prop


# NRI vs. Richness (changes from June to January)
ggplot(mapping = aes(x = richness_diff, y = NRI_diff)) + 
  geom_point(aes(color = mean_latitudes), size=4) + 
  scale_colour_viridis(name = "BCR Latitude", direction = 1) + 
  ggtitle("NRI vs. Richness Differences (June - January)") + 
  labs(x = "Species Richness Difference", y = "NRI Difference") +
  theme(text = element_text(size = 25)) +
  geom_smooth(method = "lm")

#NTI vs. Richness (changes from June to January)
ggplot(mapping = aes(x = richness_diff, y = NTI_diff)) + 
  geom_point(aes(color = mean_latitudes), size=4) + 
  scale_colour_viridis(name = "BCR Latitude", direction = 1) + 
  ggtitle("NTI vs. Richness Differences (June - January)") + 
  labs(x = "Species Richness Difference", y = "NTI Difference") +
  theme(text = element_text(size = 25)) +
  geom_smooth(method = "lm")

# NRI vs. Migratory Proportion (changes from June to January)
ggplot(mapping = aes(x = mig_prop_diff, y = NRI_diff)) + 
  geom_point(aes(color = mean_latitudes), size=4) + 
  scale_colour_viridis(name = "BCR Latitude", direction = 1) + 
  ggtitle("NRI vs. Mig. Prop. Differences (June - January)") + 
  labs(x = "Mig. Prop. Difference", y = "NRI Difference") +
  theme(text = element_text(size = 25)) +
  geom_smooth(method = "lm")

#NTI vs. Migratory Proportion (changes from June to January)
ggplot(mapping = aes(x = mig_prop_diff, y = NTI_diff)) + 
  geom_point(aes(color = mean_latitudes), size=4) + 
  scale_colour_viridis(name = "BCR Latitude", direction = 1) + 
  ggtitle("NTI vs. Mig. Prop. Differences (June - January)") + 
  theme(text = element_text(size = 25)) +
  labs(x = "Mig. Prop. Difference", y = "NTI Difference") +
  geom_smooth(method = "lm")



ggplot(filter(mean_NTI_values_visuals, variable %in% c(1:10, 23:31), BCR %in% c(12)), aes(x = richness, y = value)) + 
  geom_point(aes(color = variable)) + 
  scale_colour_viridis_d()



#showing here that there at least is significant correlation with latitude and richness. 

mean_NTI_values_visuals
test_NTI_lm <- lm(value ~ latitude + richness, data=mean_NTI_values_visuals)
summary(test_NTI_lm)

#look at correlation between migratory proportion and richness
test_NTI_lm_mig <- lm(value ~ mig_prop, data=mean_NTI_values_visuals)
summary(test_NTI_lm_mig) #significant relationship

#against all three?
test_NTI_lm_mig_2 <- lm(value ~ mig_prop + latitude + richness, data=mean_NTI_values_visuals)
summary(test_NTI_lm_mig_2) #significant relationship




#try to show something statistically here

library(lme4)
library(lmerTest)
library(sjPlot)
library(ggeffects)

#attempt doing it hierarchically

testhier <- lmer(value ~ richness + (richness | BCR), data=mean_NTI_values_visuals)
ranova(testhier)
tab_model(testhier)
as.data.frame(tab_model(testhier))



test_pred <- predict_response(testhier, terms=c("richness","BCR [sample = 33]"), type="random", interval="confidence")

#BCR-level NTI-richness association

ggplot(as.data.frame(test_pred), aes(x = x, y = predicted, group = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), color = FALSE, alpha = 0.1) +
  geom_line(color = "gray50") +
  labs(x = "Richness", y = "NTI") + 
  ggtitle("BCR-Level NTI-Richness Association") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20)) +
  geom_segment(aes(x=0, y=.779704, xend=200, yend=.779704-.004971*200), linewidth=1.1)


testhier2 <- lmer(value ~ richness + (1|variable), data=mean_NTI_values_visuals)
ranova(testhier2)
tab_model(testhier2)

testhier3 <- lmer(value ~ richness + (richness|BCR), data=mean_NRI_values_visuals)
ranova(testhier3)
tab_model(testhier3)

test_pred_2 <- predict_response(testhier3, terms=c("richness","BCR [sample = 33]"), type="random", interval="confidence")


#BCR-level NRI-richness association
#combine into one dataframe
combined_week_2_NRI <- bind_rows(week2bcrlist, .id = "runs")

head(combined_week_2_nri)

country_map + 
  geom_point(data = combined_week_2_NRI, alpha = 0.5, aes(x = x, y = y, color = ntaxa)) +
  scale_color_viridis() + 
  theme()

combined_week_26_nri <- bind_rows(week26bcrlist, .id = "runs")
country_map + 
  geom_point(data = combined_week_26_NRI, alpha = 0.5, aes(x = x, y = y, color = ntaxa)) +
  scale_color_viridis() + 
  theme()



ggplot(as.data.frame(test_pred_2), aes(x = x, y = predicted, group = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), color = FALSE, alpha = 0.1) +
  geom_line(color = "gray50") +
  labs(x = "Richness", y = "NRI") +
  ggtitle("BCR-Level NRI-Richness Association") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20)) +
  geom_segment(aes(x=0, y=-.535865, xend=200, yend=-.535865+.003055*200), linewidth=1.1)


testhier4 <- lmer(value ~ latitude + (1|variable), data=mean_NTI_values_visuals)
ranova(testhier4)
tab_model(testhier4)

testhier5 <- lmer(value ~ latitude + (1|variable), data=mean_NRI_values_visuals)
ranova(testhier5)
tab_model(testhier5)

testhier6 <- lmer(value ~ richness + mig_prop + (1 | BCR), data = mean_NTI_values_visuals)
ranova(testhier6)
tab_model(testhier6)

testhier7 <- lmer(value ~ richness + mig_prop + (1 | BCR), data = mean_NRI_values_visuals)
ranova(testhier7)
tab_model(testhier7)

testhier8 <- lmer(value ~ mig_prop + (1|BCR), data=mean_NTI_values_visuals)
ranova(testhier8)
tab_model(testhier8)

testhier9 <- lmer(value ~ mig_prop + (1|BCR), data=mean_NRI_values_visuals)
ranova(testhier9)
tab_model(testhier9)



### model comparison

library(AICcmodavg)
###NTIs
#all lmers
lmer_NTI_richness_migprop_BCR <- lmer(value ~ richness + mig_prop + (1 | BCR), data = mean_NTI_values_visuals)
lmer_NTI_richness_BCR <- lmer(value ~ richness + (1 | BCR), data = mean_NTI_values_visuals)
lmer_NTI_mig_prop_BCR <- lmer(value ~ mig_prop + (1 | BCR), data = mean_NTI_values_visuals)
lmer_NTI_mig_prop_interaction_BCR <- lmer(value ~ richness + mig_prop + richness*mig_prop + (1 | BCR), data = mean_NTI_values_visuals)
lmer_NTI_null <- lmer(value ~ (1 | BCR), data = mean_NTI_values_visuals)


NTI_models <- list(lmer_NTI_richness_migprop_BCR, lmer_NTI_mig_prop_interaction_BCR,
                   lmer_NTI_richness_BCR, lmer_NTI_mig_prop_BCR, lmer_NTI_null)

NTI_model.names <- c('lmer_NTI_richness_migprop_BCR', 'lmer_NTI_mig_prop_interaction_BCR',
                     'lmer_NTI_richness_BCR', 'lmer_NTI_mig_prop_BCR','lmer_NTI_null')

aictab(cand.set = NTI_models, modnames = NTI_model.names)
tab_model(lmer_NTI_richness_migprop_BCR)
tab_model(lmer_NTI_mig_prop_interaction_BCR)
tab_model(lmer_NTI_richness_BCR)
tab_model(lmer_NTI_mig_prop_BCR)


###NRIs
#all lmers
lmer_NRI_richness_migprop_BCR <- lmer(value ~ richness + mig_prop + (1 | BCR), data = mean_NRI_values_visuals)
lmer_NRI_richness_BCR <- lmer(value ~ richness + (1 | BCR), data = mean_NRI_values_visuals)
lmer_NRI_mig_prop_BCR <- lmer(value ~ mig_prop + (1 | BCR), data = mean_NRI_values_visuals)
lmer_NRI_mig_prop_interaction_BCR <- lmer(value ~ richness + mig_prop + richness*mig_prop + (1 | BCR), data = mean_NRI_values_visuals)
lmer_NRI_null <- lmer(value ~ (1 | BCR), data = mean_NRI_values_visuals)

NRI_models <- list(lmer_NRI_richness_migprop_BCR, lmer_NRI_mig_prop_interaction_BCR,
                   lmer_NRI_richness_BCR, lmer_NRI_mig_prop_BCR, lmer_NRI_null)

NRI_model.names <- c('lmer_NRI_richness_migprop_BCR', 'lmer_NRI_mig_prop_interaction_BCR',
                     'lmer_NRI_richness_BCR', 'lmer_NRI_mig_prop_BCR', 'lmer_NRI_null')

aictab(cand.set = NRI_models, modnames = NRI_model.names)

tab_model(lmer_NRI_richness_migprop_BCR)
tab_model(lmer_NRI_mig_prop_interaction_BCR)
tab_model(lmer_NRI_richness_BCR)
tab_model(lmer_NRI_mig_prop_BCR)


###now see what happens with just lms -- not broken into BCRs
#NTI
lm_NTI_richness_migprop <- lm(value ~ richness + mig_prop, data = mean_NTI_values_visuals)
lm_NTI_richness <- lm(value ~ richness, data = mean_NTI_values_visuals)
lm_NTI_mig_prop <- lm(value ~ mig_prop, data = mean_NTI_values_visuals)
lm_NTI_mig_prop_interaction <- lm(value ~ richness + mig_prop + richness*mig_prop, data = mean_NTI_values_visuals)

NTI_models_lm <- list(lm_NTI_richness_migprop, lm_NTI_mig_prop_interaction,
                      lm_NTI_richness, lm_NTI_mig_prop)

NTI_models_lm.names <- c('lm_NTI_richness_migprop', 'lm_NTI_mig_prop_interaction',
                         'lm_NTI_richness', 'lm_NTI_mig_prop')

aictab(cand.set = NTI_models_lm, modnames = NTI_models_lm.names)

#NRI
lm_NRI_richness_migprop <- lm(value ~ richness + mig_prop, data = mean_NRI_values_visuals)
lm_NRI_richness <- lm(value ~ richness, data = mean_NRI_values_visuals)
lm_NRI_mig_prop <- lm(value ~ mig_prop, data = mean_NRI_values_visuals)
lm_NRI_mig_prop_interaction <- lm(value ~ richness + mig_prop + richness*mig_prop, data = mean_NRI_values_visuals)

NRI_models_lm <- list(lm_NRI_richness_migprop, lm_NRI_mig_prop_interaction,
                      lm_NRI_richness, lm_NRI_mig_prop)

NRI_models_lm.names <- c('lm_NRI_richness_migprop', 'lm_NRI_mig_prop_interaction',
                         'lm_NRI_richness', 'lm_NRI_mig_prop')

aictab(cand.set = NRI_models_lm, modnames = NRI_models_lm.names)




#create mixed linear effects model for NTI on BCR-level data
testmlm <- lmer(value ~ richness + (1 | BCR), data = mean_NTI_values_visuals)
summary(testmlm)
confint(testmlm)

#create mixed linear effects model for richness and mig_prop on BCR-level data

testmlm2 <- lmer(value ~ richness + mig_prop + (1 | BCR), data = mean_NTI_values_visuals)
summary(testmlm2)
testmlm3 <- lmer(value ~ richness + mig_prop + (1 | BCR), data = mean_NRI_values_visuals)
summary(testmlm3)



#estimates of random effects
coef(testmlm)$BCR %>% head(5)

library(merTools)
library(lattice)

#run and plot estimated linear effects
REsim(testmlm) 
plotREsim(REsim(testmlm))
?plotREsim()
?lmer()

#dotplot of random effects (there's geographic signal)
dotplot(ranef(testmlm))
random_effect_intercepts <- as.data.frame(ranef(testmlm))


#start a new dataframe that makes BCR data modifiable
allbcrs_extra_data <- allbcrs
allbcrs_extra_data$random_intercepts <- rep(NA, nrow(allbcrs_extra_data))
allbcrs_extra_data[allbcrs_extra_data$BCR == 5,"random_intercepts"]

for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"random_intercepts"] <- random_effect_intercepts[random_effect_intercepts$grp == i,"condval"]
}


ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = random_intercepts)) + 
  coord_sf(xlim = c(-179, -55),
           ylim = c(25, 80)) + 
  scale_fill_gradient2("Random Intercept",midpoint=0, low="red", mid="white",
                       high="blue")

###interpretation: blue = lower NTI in summer = more summer overdispersion, more winter clustering
#                 red = lower NTI in winter = more winter overdispersion, more summer clustering
#                 migrants bringing their overdispersion northwards in the summer




#replicate with NRI
#create mixed linear effects model for NRI on BCR-level data
NRI_mlm <- lmer(value ~ richness + (1 | BCR), data = mean_NRI_values_visuals)
?lmer()
dotplot(ranef(NRI_mlm))
plotREsim(REsim(NRI_mlm))
random_effect_intercepts_NRI <- as.data.frame(ranef(NRI_mlm))

#add to BCR data
allbcrs_extra_data$random_intercepts_NRI <- rep(NA, nrow(allbcrs_extra_data))
for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"random_intercepts_NRI"] <- random_effect_intercepts_NRI[random_effect_intercepts_NRI$grp == i,"condval"]
}

ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = random_intercepts_NRI)) + 
  coord_sf(xlim = c(-179, -55),
           ylim = c(25, 80)) + 
  scale_fill_gradient2("Random Intercept",midpoint=0, low="red", mid="white",
                       high="blue")


###interpretation: blue = lower NTI in summer = more summer overdispersion, more winter clustering
#                 red = lower NTI in winter = more winter overdispersion, more summer clustering
#                 a cloudier picture overall, but some element of the pattern being reversed


# plot richness diff, NRI diff, and NTI diff on the same plot space
allbcrs_extra_data$NRI_diff <- rep(NA, nrow(allbcrs_extra_data))
allbcrs_extra_data$NTI_diff <- rep(NA, nrow(allbcrs_extra_data))
allbcrs_extra_data$richness_diff <- rep(NA, nrow(allbcrs_extra_data))

for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"NRI_diff"] <- NRI_diff[i-4]
}
for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"NTI_diff"] <- NTI_diff[i-4]
}
for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"richness_diff"] <- richness_diff[i-4]
}

#### manuscript maps here

#plot NRI difference
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = NRI_diff)) + 
  coord_sf(xlim = c(-179, -55),
           ylim = c(25, 80)) + 
  scale_fill_gradient2("NRI Difference",midpoint=0, low="red", mid="white",
                       high="blue")

#plot NTI difference
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = NTI_diff)) + 
  coord_sf(xlim = c(-179, -55),
           ylim = c(25, 80)) + 
  scale_fill_gradient2("NTI Difference",midpoint=0, low="red", mid="white",
                       high="blue")

#plot richness difference
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = richness_diff)) + 
  coord_sf(xlim = c(-179, -55),
           ylim = c(25, 80)) + 
  scale_fill_gradient2("Richness Difference",midpoint=0, low="red", mid="white",
                       high="blue")

#what can I do with random slopes?

#adding random slope to NTI data
NTI_mlm_random_slope <- lmer(value ~ richness + (richness|BCR), data = mean_NTI_values_visuals)
dotplot(ranef(NTI_mlm_random_slope),
        scales = list(x =list(relation = 'free')))

random_effect_slope_intercepts_NTI <- as.data.frame(ranef(NTI_mlm_random_slope))
random_effect_slope_intercepts_NTI

random_effect_slope_intercepts_NTI[random_effect_slope_intercepts_NTI$grp == i & 
                                     random_effect_slope_intercepts_NTI$term == "richness","condval"]

allbcrs_extra_data$NTI_slope <- rep(NA, nrow(allbcrs_extra_data))
for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"NTI_slope"] <- random_effect_slope_intercepts_NTI[random_effect_slope_intercepts_NTI$grp == i & 
                                                                                                      random_effect_slope_intercepts_NTI$term == "richness","condval"]
}

#plotting random slopes NTI data
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = NTI_slope)) + 
  coord_sf(xlim = c(-179, -55),
           ylim = c(25, 80)) + 
  scale_fill_gradient2("Random Slope",midpoint=0, low="red", mid="white",
                       high="blue")


#adding random slope to NRI data
NRI_mlm_random_slope <- lmer(value ~ richness + (richness|BCR), data = mean_NRI_values_visuals)
dotplot(ranef(NRI_mlm_random_slope),
        scales = list(x =list(relation = 'free')))

random_effect_slope_intercepts_NRI <- as.data.frame(ranef(NRI_mlm_random_slope))
random_effect_slope_intercepts_NRI

allbcrs_extra_data$NRI_slope <- rep(NA, nrow(allbcrs_extra_data))
for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"NRI_slope"] <- random_effect_slope_intercepts_NRI[random_effect_slope_intercepts_NRI$grp == i & 
                                                                                                      random_effect_slope_intercepts_NRI$term == "richness","condval"]
}

#plotting random slopes NRI data
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = NRI_slope)) + 
  coord_sf(xlim = c(-179, -55),
           ylim = c(25, 80)) + 
  scale_fill_gradient2("Random Slope",midpoint=0, low="red", mid="white",
                       high="blue")


#and attempt visualizing the mixed effect regressions
ggplot(mean_NTI_values_visuals, aes(x = richness, y = value)) + 
  geom_point(aes(colour = richness)) + 
  facet_wrap(~BCR) + 
  scale_colour_viridis(direction = 1) + 
  ggtitle("Annual NTI Distribution by BCR") + 
  labs(x = "Week", y = "NTI")

random_effect_slope_intercepts_NTI


### see what happens if you omit migration
NTI_mlm_random_slope_nomig <- lmer(value ~ richness + (richness|BCR), data = mean_NTI_values_visuals[mean_NTI_values_visuals$variable %in% c(1:10, 21:33, 46:52),])
dotplot(ranef(NTI_mlm_random_slope_nomig),
        scales = list(x =list(relation = 'free')))

random_effect_slope_intercepts_NTI_nomig <- as.data.frame(coef(NTI_mlm_random_slope_nomig))
random_effect_slope_intercepts_NTI_nomig



#(and with NRI)

NRI_mlm_random_slope_nomig <- lmer(value ~ richness + (richness|BCR), data = mean_NRI_values_visuals[mean_NRI_values_visuals$variable %in% c(1:10, 21:33, 46:52),])
dotplot(ranef(NRI_mlm_random_slope_nomig),
        scales = list(x =list(relation = 'free')))

random_effect_slope_intercepts_NRI_nomig <- as.data.frame(ranef(NRI_mlm_random_slope_nomig))

allbcrs_extra_data$NRI_slope_nomig <- rep(NA, nrow(allbcrs_extra_data))
for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"NRI_slope_nomig"] <- random_effect_slope_intercepts_NRI_nomig[random_effect_slope_intercepts_NRI_nomig$grp == i & 
                                                                                                                  random_effect_slope_intercepts_NRI_nomig$term == "richness","condval"]
}





#trying to print everything

library(ggeffects)
library(sjPlot)
library(tibble)

coef(NTI_mlm_random_slope_nomig)


# see group coefficients
model_coefs <- coef(NTI_mlm_random_slope_nomig)$BCR %>% 
  rename(Intercept = `(Intercept)`, Slope = richness) %>% 
  rownames_to_column("BCR")

# see coefficients
model_coefs

#left join
mean_NTI_values_visuals$BCR <- as.character(mean_NTI_values_visuals$BCR)
mean_NTI_values_visuals_nomig <- mean_NTI_values_visuals[mean_NTI_values_visuals$variable %in% c(1:10, 21:33, 46:52),]
mean_NTI_values_visuals_nomig_rans <- left_join(mean_NTI_values_visuals_nomig, model_coefs, by = "BCR")


#### plot individual random effect slope and intercept values
model_coef_plot <- ggplot(data = mean_NTI_values_visuals_nomig_rans, 
                          mapping = aes(x = richness, 
                                        y = value, 
                                        colour = BCR)
) +
  geom_point(na.rm = T, alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, 
                  slope = Slope,
                  colour = BCR
  ),
  linewidth = 0.75
  ) +
  facet_wrap(~BCR)

model_coef_plot %+% mean_NTI_values_visuals_nomig_rans

#fill this data into all BCR data
allbcrs_extra_data$NTI_slope_nomig <- rep(NA, nrow(allbcrs_extra_data))
for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"NTI_slope_nomig"] <- model_coefs[model_coefs$BCR == i, "Slope"]
}


#plot
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = NTI_slope_nomig)) + 
  coord_sf(xlim = c(-179, -55),
           ylim = c(25, 80)) + 
  scale_fill_gradient2("Random Slope",midpoint=0, low="red", mid="white",
                       high="blue")


ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = richness_diff)) + 
  coord_sf(xlim = c(-179, -55),
           ylim = c(25, 80)) + 
  scale_fill_gradient2("Richness Diff",midpoint=0, low="red", mid="white",
                       high="blue")


cor(x=allbcrs_extra_data$NTI_slope_nomig, y=allbcrs_extra_data$richness_diff, use="complete.obs")

summary(lm(model_coefs$Slope ~ richness_diff))

richness_diff
model_coefs$Slope

plot(richness_diff ~ NTI_slope_nomig, data=allbcrs_extra_data)



##### replicate with continental-scale regional species pool

#import each NRI and NTI table
for(i in 1:52){
  assign(paste0("week_", i, "_NRI"), read.csv(paste0("output/tables/", i, "/passerine.ses.mpd.", i, ".csv"), row.names=1))
  assign(paste0("week_", i, "_NTI"), read.csv(paste0("output/tables/", i, "/passerine.ses.mntd.", i, ".csv"), row.names=1))
}

head(week_1_NRI)
head(combined_week_26_NRI)

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

#use agg spatraster of minimum length df
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


#apply BCR values to NRI list (warning! takes long processing time)
#for(i in 1:length(NRIlist)){
#  for(j in 1:nrow(NRIlist[[i]])){
#    if((length(allbcrrows[allbcrrows$rows == rownames((NRIlist[[i]])[j,]),"BCR"]) == 1) == TRUE){
#      NRIlist[[i]][j,12] <- allbcrrows[allbcrrows$rows == rownames((NRIlist[[i]])[j,]),"BCR"]
#    } else {
#      print("Pass")
#    }
#  }
#}

#make one dataframe
#NRIlistdf <- do.call(rbind, unname(NRIlist))

#write.csv(NRIlistdf, "output/tables/all_NRI_data.csv")
NRIlistdf <- read.csv("output/tables/all_NRI_data.csv")


## now go on with continental NTI
#create a list of all NRI dataframes 
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

#apply BCR values to NTI list (warning! takes long processing time)
#for(i in 1:length(NTIlist)){
# for(j in 1:nrow(NTIlist[[i]])){
#   if((length(allbcrrows[allbcrrows$rows == rownames((NTIlist[[i]])[j,]),"BCR"]) == 1) == TRUE){
#     NTIlist[[i]][j,12] <- allbcrrows[allbcrrows$rows == rownames((NTIlist[[i]])[j,]),"BCR"]
#   } else {
#     print("Pass")
#   }
# }
#}

#make one dataframe
#NTIlistdf <- do.call(rbind, unname(NTIlist))

#write.csv(NTIlistdf, "output/tables/all_NTI_data.csv")
NTIlistdf <- read.csv("output/tables/all_NTI_data.csv")
head(NTIlistdf)


#### add latitude and longitude data to NTIlistdf and NRIlistdf
NTIlistdf.rows <- as.numeric(NTIlistdf$X)
NRIlistdf.rows <- as.numeric(NRIlistdf$X)


NTIlistdf$x <- rep(NA, nrow(NTIlistdf))
NTIlistdf$y <- rep(NA, nrow(NTIlistdf))

NRIlistdf$x <- rep(NA, nrow(NRIlistdf))
NRIlistdf$y <- rep(NA, nrow(NRIlistdf))


for(i in 1:nrow(NRIlistdf)){
  NRIlistdf$x[i] <- xyFromCell(eawpew_26_agg, NRIlistdf.rows[i])[,1]
  NRIlistdf$y[i] <- xyFromCell(eawpew_26_agg, NRIlistdf.rows[i])[,2]
}

for(i in 1:nrow(NTIlistdf)){
  NTIlistdf$x[i] <- xyFromCell(eawpew_26_agg, NTIlistdf.rows[i])[,1]
  NTIlistdf$y[i] <- xyFromCell(eawpew_26_agg, NTIlistdf.rows[i])[,2]
}


NRIlistdf$latitude <- sf_project(from = crs, to = st_crs(4326), NRIlistdf[,c("x","y")])[,2]
NRIlistdf$longitude <- sf_project(from = crs, to = st_crs(4326), NRIlistdf[,c("x","y")])[,1]

NTIlistdf$latitude <- sf_project(from = crs, to = st_crs(4326), NTIlistdf[,c("x","y")])[,2]
NTIlistdf$longitude <- sf_project(from = crs, to = st_crs(4326), NTIlistdf[,c("x","y")])[,1]

#### add migratory proportion
NRIlistdf$proportion <-  rep(NA, nrow(NRIlistdf))
NTIlistdf$proportion <-  rep(NA, nrow(NTIlistdf))

#load all_mig_props from community_migratory_proportion.R
#will take a long time
#occasionally you hit an NA for location and have to skip it
#load week command from below
for(i in 1:nrow(NRIlistdf)){ #changed value based on where it left off
  NRIlistdf[i,"proportion"] <- all_mig_props[[NRIlistdf[i,"week"]]][all_mig_props[[NRIlistdf[i,"week"]]]$community == NRIlistdf[i,"X"],"proportion"]
}


### populate same structure

#NRI continental
#create table of mean NRI values
BCRs <- c(5:37)
weeks <- c(1:52)
mean_NRI_values_matrix_continental <- matrix(NA, ncol = length(weeks), nrow = length(BCRs))
col_names <- as.character(weeks)

# Convert the matrix to a data frame with column names
mean_NRI_values_continental <- data.frame(matrix = mean_NRI_values_matrix_continental)
colnames(mean_NRI_values_continental) <- col_names
rownames(mean_NRI_values_continental) <- BCRs
head(mean_NRI_values_continental)

#add weeks to the NRI list df
NRIlistdf$week <- rep(NA, nrow(NRIlistdf))
NRIlistdf[c(1:25000),"week"] <- 1
for(i in 1:51){
  NRIlistdf[(25000*i)+c(1:25000),"week"] <- i+1
}

mean_NRI_values_continental$BCR <- BCRs

#populate the mean NRI values table
for(i in BCRs){ 
  for(j in weeks){
    mean_NRI_values_continental[i-4, j] <- mean(na.omit(NRIlistdf[NRIlistdf$week == j & NRIlistdf$BCR == i, "NRI"]))
  }
}

head(mean_NRI_values_continental)
mean_NRI_values_continental$latitude <- mean_latitudes

#NTI continental
BCRs <- c(5:37)
weeks <- c(1:52)
mean_NTI_values_matrix_continental <- matrix(NA, ncol = length(weeks), nrow = length(BCRs))
col_names <- as.character(weeks)

# Convert the matrix to a data frame with column names
mean_NTI_values_continental <- data.frame(matrix = mean_NTI_values_matrix_continental)
colnames(mean_NTI_values_continental) <- col_names
rownames(mean_NTI_values_continental) <- BCRs
head(mean_NTI_values_continental)

#add weeks to the NTI list df
NTIlistdf$week <- rep(NA, nrow(NTIlistdf))
NTIlistdf[c(1:25000),"week"] <- 1
for(i in 1:51){
  NTIlistdf[(25000*i)+c(1:25000),"week"] <- i+1
}

mean_NTI_values_continental$BCR <- BCRs

#populate the mean NTI values table
for(i in BCRs){ 
  for(j in weeks){
    mean_NTI_values_continental[i-4, j] <- mean(na.omit(NTIlistdf[NTIlistdf$week == j & NTIlistdf$BCR == i, "NTI"]))
  }
}

head(mean_NTI_values_continental)

mean_NTI_values_continental$latitude <- mean_latitudes

#should be able to use the same values for species diversity but can go back to redo if need be

#NRI -- continental
mean_NRI_values_visuals_continental <- reshape2::melt(mean_NRI_values_continental, id.vars=c("BCR", "latitude"))
mean_NRI_values_visuals_continental$richness <- weekly_bcr_diversity_values
mean_NRI_values_visuals_continental$mig_prop <- weekly_bcr_mig_prop_values

#NTI -- continental
mean_NTI_values_visuals_continental <- reshape2::melt(mean_NTI_values_continental, id.vars=c("BCR", "latitude"))
mean_NTI_values_visuals_continental$richness <- weekly_bcr_diversity_values
mean_NTI_values_visuals_continental$mig_prop <- weekly_bcr_mig_prop_values

#create NRI, NTI, and richness difference vectors
week1_NRIs_continental <- mean_NRI_values_visuals_continental[mean_NRI_values_visuals_continental$variable == 1,"value"]
week1_NTIs_continental <- mean_NTI_values_visuals_continental[mean_NTI_values_visuals_continental$variable == 1,"value"]
week1_richness_continental <- mean_NRI_values_visuals_continental[mean_NRI_values_visuals_continental$variable == 1,"richness"]
week26_NRIs_continental <- mean_NRI_values_visuals_continental[mean_NRI_values_visuals_continental$variable == 26,"value"]
week26_NTIs_continental <- mean_NTI_values_visuals_continental[mean_NTI_values_visuals_continental$variable == 26,"value"]
week26_richness_continental <- mean_NRI_values_visuals_continental[mean_NRI_values_visuals_continental$variable == 26,"richness"]
NRI_diff_continental <- week26_NRIs_continental - week1_NRIs_continental
NTI_diff_continental <- week26_NTIs_continental - week1_NTIs_continental
richness_diff_continental <- week26_richness_continental - week1_richness_continental
week1_mig_prop_continental <- mean_NRI_values_visuals_continental[mean_NRI_values_visuals_continental$variable == 1,"mig_prop"]
week26_mig_prop_continental <- mean_NRI_values_visuals_continental[mean_NRI_values_visuals_continental$variable == 26,"mig_prop"]
mig_prop_diff_continental <- week26_mig_prop_continental - week1_mig_prop_continental


# NRI vs. Richness (changes from June to January)
ggplot(mapping = aes(x = richness_diff_continental, y = NRI_diff_continental)) + 
  geom_point(aes(color = mean_latitudes), size=4) + 
  scale_colour_viridis(name = "BCR Latitude", direction = 1) + 
  ggtitle("NRI vs. Richness Differences (June - January)") + 
  labs(x = "Species Richness Difference", y = "NRI Difference") +
  theme(text = element_text(size = 25)) +
  geom_smooth(method = "lm")

#NTI vs. Richness (changes from June to January)
ggplot(mapping = aes(x = richness_diff_continental, y = NTI_diff_continental)) + 
  geom_point(aes(color = mean_latitudes), size=4) + 
  scale_colour_viridis(name = "BCR Latitude", direction = 1) + 
  ggtitle("NTI vs. Richness Differences (June - January)") + 
  labs(x = "Species Richness Difference", y = "NTI Difference") +
  theme(text = element_text(size = 25)) +
  geom_smooth(method = "lm")

# NRI vs. Migratory Proportion (changes from June to January)
ggplot(mapping = aes(x = mig_prop_diff_continental, y = NRI_diff_continental)) + 
  geom_point(aes(color = mean_latitudes), size=4) + 
  scale_colour_viridis(name = "BCR Latitude", direction = 1) + 
  ggtitle("NRI vs. Mig. Prop. Differences (June - January)") + 
  labs(x = "Mig. Prop. Difference", y = "NRI Difference") +
  theme(text = element_text(size = 25)) +
  geom_smooth(method = "lm")

#NTI vs. Migratory Proportion (changes from June to January)
ggplot(mapping = aes(x = mig_prop_diff_continental, y = NTI_diff_continental)) + 
  geom_point(aes(color = mean_latitudes), size=4) + 
  scale_colour_viridis(name = "BCR Latitude", direction = 1) + 
  ggtitle("NTI vs. Mig. Prop. Differences (June - January)") + 
  labs(x = "Mig. Prop. Difference", y = "NTI Difference") +
  theme(text = element_text(size = 25)) +
  geom_smooth(method = "lm")

### model comparison

###NTIs
#all lmers
lmer_NTI_richness_migprop_BCR_continental <- lmer(value ~ richness + mig_prop + (1 | BCR), data = mean_NTI_values_visuals_continental)
lmer_NTI_richness_BCR_continental <- lmer(value ~ richness + (1 | BCR), data = mean_NTI_values_visuals_continental)
lmer_NTI_mig_prop_BCR_continental <- lmer(value ~ mig_prop + (1 | BCR), data = mean_NTI_values_visuals_continental)
lmer_NTI_mig_prop_interaction_BCR_continental <- lmer(value ~ richness + mig_prop + richness*mig_prop + (1 | BCR), data = mean_NTI_values_visuals_continental)
lmer_NTI_null_continental <- lmer(value ~ (1 | BCR), data = mean_NTI_values_visuals_continental)


NTI_models_continental <- list(lmer_NTI_richness_migprop_BCR_continental, lmer_NTI_mig_prop_interaction_BCR_continental,
                               lmer_NTI_richness_BCR_continental, lmer_NTI_mig_prop_BCR_continental, lmer_NTI_null_continental)

NTI_model.names_continental <- c('lmer_NTI_richness_migprop_BCR_continental', 'lmer_NTI_mig_prop_interaction_BCR_continental',
                                 'lmer_NTI_richness_BCR_continental', 'lmer_NTI_mig_prop_BCR_continental','lmer_NTI_null_continental')

aictab(cand.set = NTI_models_continental, modnames = NTI_model.names_continental)
tab_model(lmer_NTI_mig_prop_BCR_continental)
tab_model(lmer_NTI_richness_BCR_continental)
tab_model(lmer_NTI_richness_migprop_BCR_continental)
tab_model(lmer_NTI_mig_prop_interaction_BCR_continental)


###NRIs
#all lmers
lmer_NRI_richness_migprop_BCR_continental <- lmer(value ~ richness + mig_prop + (1 | BCR), data = mean_NRI_values_visuals_continental)
lmer_NRI_richness_BCR_continental <- lmer(value ~ richness + (1 | BCR), data = mean_NRI_values_visuals_continental)
lmer_NRI_mig_prop_BCR_continental <- lmer(value ~ mig_prop + (1 | BCR), data = mean_NRI_values_visuals_continental)
lmer_NRI_mig_prop_interaction_BCR_continental <- lmer(value ~ richness + mig_prop + richness*mig_prop + (1 | BCR), data = mean_NRI_values_visuals_continental)
lmer_NRI_null_continental <- lmer(value ~ (1 | BCR), data = mean_NRI_values_visuals_continental)

NRI_models_continental <- list(lmer_NRI_richness_migprop_BCR_continental, lmer_NRI_mig_prop_interaction_BCR_continental,
                               lmer_NRI_richness_BCR_continental, lmer_NRI_mig_prop_BCR_continental, lmer_NRI_null_continental)

NRI_model.names_continental <- c('lmer_NRI_richness_migprop_BCR_continental', 'lmer_NRI_mig_prop_interaction_BCR_continental',
                                 'lmer_NRI_richness_BCR_continental', 'lmer_NRI_mig_prop_BCR_continental', 'lmer_NRI_null_continental')

aictab(cand.set = NRI_models_continental, modnames = NRI_model.names_continental)

tab_model(lmer_NRI_mig_prop_interaction_BCR_continental)
tab_model(lmer_NRI_mig_prop_BCR_continental)
tab_model(lmer_NRI_richness_BCR_continental)
tab_model(lmer_NRI_richness_migprop_BCR_continental)

### correlation of richness and mig prop
head(mean_NRI_values_visuals_continental)
plot(richness ~ mig_prop, data=mean_NRI_values_visuals_continental[mean_NRI_values_visuals_continental$mig_prop > 0,])
mean_NRI_values_visuals_continental[mean_NRI_values_visuals_continental$mig_prop == 0,]


tab_model(lm(richness ~ mig_prop, data=mean_NRI_values_visuals_continental[mean_NRI_values_visuals_continental$mig_prop > 0,]))
cor.test(mean_NRI_values_visuals_continental[mean_NRI_values_visuals_continental$mig_prop > 0,"richness"],
         mean_NRI_values_visuals_continental[mean_NRI_values_visuals_continental$mig_prop > 0,"mig_prop"])

weekly_bcr_mig_prop_values[weekly_bcr_mig_prop_values > 0]



### try to recreate the seasonal continental differences maps
allbcrs_extra_data$NRI_diff_continental <- rep(NA, nrow(allbcrs_extra_data))
allbcrs_extra_data$NTI_diff_continental <- rep(NA, nrow(allbcrs_extra_data))
allbcrs_extra_data$mig_prop_diff <- rep(NA, nrow(allbcrs_extra_data))

for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"NRI_diff_continental"] <- NRI_diff_continental[i-4]
}
for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"NTI_diff_continental"] <- NTI_diff_continental[i-4]
}
for(i in 5:37){
  allbcrs_extra_data[allbcrs_extra_data$BCR == i,"mig_prop_diff"] <- mig_prop_diff_continental[i-4]
}


#### manuscript maps here
library(mapproj)

crs_102004 <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

allbcrs_extra_data <- st_transform(allbcrs_extra_data, crs = crs_102004)
allbcrs_extra_data_coords <- st_bbox(allbcrs_extra_data)

#plot NRI difference
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = NRI_diff)) + 
  coord_sf(crs = crs_102004, 
           xlim = c(-3550000, 3000000), 
           ylim = c(-1700000, 4000000)) + 
  ggtitle("Seasonal NRI Difference (June - January), Regional Pool") +
  theme(text = element_text(size = 20)) +
  scale_fill_gradient2("NRI Diff.",midpoint=0, low="red", mid="white",
                       high="blue")

#plot NTI difference
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = NTI_diff)) + 
  coord_sf(crs = crs_102004, 
           xlim = c(-3550000, 3000000), 
           ylim = c(-1700000, 4000000)) + 
  ggtitle("Seasonal NTI Difference (June - January), Regional Pool") +
  theme(text = element_text(size = 20)) +
  scale_fill_gradient2("NTI Diff.",midpoint=0, low="red", mid="white",
                       high="blue")

#plot NRI difference (continental)
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = NRI_diff_continental)) + 
  coord_sf(crs = crs_102004, 
           xlim = c(-3550000, 3000000), 
           ylim = c(-1700000, 4000000)) + 
  ggtitle("Seasonal NRI Difference (June - January), Continental Pool") +
  theme(text = element_text(size = 20)) +
  scale_fill_gradient2("NRI Diff.",midpoint=0, low="red", mid="white",
                       high="blue")

#plot NTI difference (continental)
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = NTI_diff_continental)) + 
  coord_sf(crs = crs_102004, 
           xlim = c(-3550000, 3000000), 
           ylim = c(-1700000, 4000000)) + 
  ggtitle("Seasonal NTI Difference (June - January), Continental Pool") +
  theme(text = element_text(size = 20)) +
  scale_fill_gradient2("NTI Diff.",midpoint=0, low="red", mid="white",
                       high="blue")


#plot richness difference
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = richness_diff)) + 
  coord_sf(crs = crs_102004, 
           xlim = c(-3550000, 3000000), 
           ylim = c(-1700000, 4000000)) + 
  ggtitle("Seasonal Species Richness Difference (June - January)") +
  theme(text = element_text(size = 20)) +
  scale_fill_gradient2("Richness Diff.",midpoint=0, low="red", mid="white",
                       high="blue")

#plot migratory proportion difference
ggplot() + 
  geom_sf(data = allbcrs_extra_data, aes(fill = mig_prop_diff)) + 
  coord_sf(crs = crs_102004, 
           xlim = c(-3550000, 3000000), 
           ylim = c(-1700000, 4000000)) + 
  ggtitle("Seasonal Migratory Proportion Difference (June - January)") +
  theme(text = element_text(size = 20)) +
  scale_fill_gradient2("Mig. Prop. Diff.",midpoint=0, low="red", mid="white",
                       high="blue")


