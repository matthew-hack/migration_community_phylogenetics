# 1. SDM download and community formation
# This script downloads the eBird Status and Trends data for each species included in the dataset.
# The eBirdst data is partitioned into weekly distributional rasters.
# These weekly species rasters are aggregated into weekly bird community matrices.

#load required libraries
required.libraries <- c("ebirdst","terra","sf","dplyr","tidyr","rnaturalearth",
                        "geodata","ggplot2","fields","stringr","TreeTools",
                        "viridis","rgdal","raster","devtools","data.table",
                        "ape","picante","pez","phytools","vegan","adephylo",
                        "phylobase","geiger","mvMORPH","OUwie","hisse",
                        "BAMMtools","phylosignal","Biostrings","kableExtra",
                        "betapart","gridExtra","reshape2")

lapply(required.libraries, require, character.only = TRUE)

extract <- terra::extract
project <- terra::project

# load the mapping parameters to generate necessary CRS
path <- get_species_path("example_data")
fac_parameters <- load_fac_map_parameters(path)
crs <- fac_parameters$custom_projection

##### DOWNLOAD eBIRD STATUS AND TRENDS DATA #####

# This data was downloaded using the 2021 eBirdst packages; a newer version should be substituted for newer downloads.

#read in CSV defining species codes and scientific names
species_codes <- read.csv("species_codes.csv")

#define species whose data is to be downloaded from eBirdst
species_names <- species_codes$ebirdst_code

#create vector of scientific names of the above species for interfacing with phylogenetic trees downstream
species_passerines <- species_codes$jetz_tree_latin_name

#download complete eBirdst data package for each species
for(i in 1:length(species_names)){
 str_sub(species_names[i]) <- ebirdst_download(species = species_names[i]) #requires a valid access key set up on eBird Status and Trends
}
#this download command generates a folder that can be located using the command ebirdst_data_dir(). File can be moved to path "data/ebirdst/2021"
#each species' complete status & trends data is downloaded into a separate folder with the name given in the species_names raster

#load example raster
ovenbi1_abd_seasonal_mean <- load_raster("data/ebirdst/2021/ovenbi1", product = "abundance", 
                                         period = "seasonal", metric = "mean", 
                                         resolution = "hr")

#use example raster to crop to US boundaries
us <- ne_countries(country = "united states of america", type = "countries", returnclass = "sf") %>% 
  st_transform(st_crs(ovenbi1_abd_seasonal_mean))

##### PROCESS WEEKLY SPECIES DISTRIBUTIONAL RASTERS #####

#run command that generates weekly aggregated csvs representing presence and absence of each species in the dataset
### WARNING: takes a very long time to run on initial pass (as long as 2 days depending on computer's processing power)
for(j in 1:52){ #iterates over each week of the year
 for(i in 1:length(species_names)){ #iterates across each species
   if(file.exists(paste0("output/species_data/", j, "/", species_names[i], "_", j, "_agg.csv")) == FALSE){ #initializes processing loop if aggregated weekly species raster does not yet exist
     assign(paste0(species_names[i], "_", j), subset(load_raster(paste0("data/ebirdst/2021/", species_names[i]), product = "abundance", resolution = "hr"), j)) #load abundance raster during iterated week of the year
     assign(paste0(species_names[i], "_", j, "_us"), crop(get(paste0(species_names[i], "_", j)), us)) #crop rasters
     assign(paste0(species_names[i], "_", j, "_proj"), project(get(paste0(species_names[i], "_", j, "_us")), crs, method="near")) #reproject rasters to CRS used by eBirdst 
     assign(paste0(species_names[i], "_", j, "_agg"), aggregate(get(paste0(species_names[i], "_", j, "_proj")), fact=2)) #aggregate rasters into appropriate scale for community analyses.
     assign(paste0(species_names[i], "_", j, "_agg_df"), as.data.frame(get(paste0(species_names[i], "_", j, "_agg")))) #convert rasters into dataframe format so they can be exported as CSVs
     write.csv(get(paste0(species_names[i], "_", j, "_agg_df")), file = paste0("output/species_data/", j, "/", species_names[i], "_", j, "_agg.csv")) #export aggregated weekly species rasters as CSVs
   }
   else { #skips regeneration of csv if file already exists
     print("File exists")
   }
 }
}

##### TRIM SONGBIRD DATASETS #####

#Each weekly eBirdst output varies slightly in scope; run these commands to trim larger weekly species csvs to only include communities for which data exists for each species.
#Rerun each week separately. Below represents an example week (week 45), which should be adjusted for the output for each week of the year.
#Warning: process can take at least an hour for each week.

#create list consisting of all species datasets in specified week
all_week_45_datasets <- rep(NA, length(species_names)) 
#initialize vector representing number of communities encompassed by each species' weekly eBirdst coverage
week_45_number_of_rows <- rep(NA, length(species_names)) 
#run loop to extract number of communities in each dataset
for(i in 286:length(species_names)){  #initialize loop to generate list of all species datasets in specified week, and extract the number of communities represented in each dataset
  assign(paste0(species_names[i], "_45_agg_df"), read.csv(paste0("output/species_data/45/", species_names[i], "_45_agg.csv"), row.names = 1)) #load all species datasets for given week
  all_week_45_datasets[i] <- get(paste0(species_names[i], "_45_agg_df")) #add species' datasets to list
  week_45_number_of_rows[i] <- nrow(as.data.frame(all_week_45_datasets[i])) #extract number of rows in each dataset and store in list
}
#identify a species (the first taxonomically) that contains the minimum number of communities encompassed by each species to trim datasets of species with greater coverage
species_names[which(week_45_number_of_rows == min(week_45_number_of_rows))[1]] #returns yellow-bellied flycatcher (yebfly) for week 45


#generate reduced distributional datasets for each species encompassing the communities for which data exists for each species
for(i in 1:length(species_names)){
  assign(paste0(species_names[i], "_45_reduced"), #create reduced distributional datasets
         get(paste0(species_names[i], "_45_agg_df"))[rownames(get(paste0(species_names[i], "_45_agg_df"))) %in% rownames(yebfly_45_agg_df),]) #trim to communities represented in each dataset; substitute individual dataset depending on results of previous step
  write.csv(get(paste0(species_names[i], "_45_reduced")), file = paste0("output/species_data/45/", species_names[i], "_45_reduced.csv")) #export reduced datasets
}

##### CREATE WEEKLY COMMUNITY MATRICES #####

#Each community spreadsheet is created individually, due to different species being used to trim the datasets. Example continued for week 45.
#Warning: process can take multiple hours per community.

#create week 45 community matrix
#load first taxonomic species at minimum community size
yebfly_45_agg_df <- read.csv("output/species_data/45/yebfly_45_agg.csv", row.names=1) #reread individual dataset if not already loaded
#create spreadsheet with number of rows = number of communities, and number of columns = number of species
community45 <- data.frame(matrix(nrow=nrow(yebfly_45_agg_df),ncol=length(species_passerines))) 
#label each row as an individually identified community
rownames(community45) <- rownames(yebfly_45_agg_df) 
#label each column as an individual species (using the scientific names, rather than the six-letter codes, will be important for downstream phylogenetic community analyses)
colnames(community45) <- species_passerines 

#populate week 45 community matrix
for(i in 1:ncol(community45)){
  community45[,i] <- get(paste0(species_names[i], "_45_reduced")) #populate each column with individual species data (species datasets created or loaded in previous step)
}
#export weekly community data matrix
data.table::fwrite(community45, file = "output/community_data/community_45.csv")
