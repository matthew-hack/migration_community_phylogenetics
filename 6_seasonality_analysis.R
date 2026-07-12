#load required libraries
required.libraries <- c("ebirdst","terra","sf","dplyr","rnaturalearth",
                        "ape","picante","vegan","ecospat","EcoSimR","betapart",
                        "ggpubr","phytools","pez","adephylo","phylobase","geiger",
                        "mvMORPH","OUwie","hisse","BAMMtools","phylosignal",
                        "Biostrings","reshape2","gridExtra","kableExtra",
                        "ggplot2","viridis","data.table","lme4","lmerTest",
                        "sjPlot","geodata","nlme","ggExtra","mgcv","sdmTMB",
                        "spdep","mapproj","ggspatial","rnaturalearthdata")

lapply(required.libraries, require, character.only = TRUE)

project <- terra::project 

#read in data
NRI_turnover_ratios <- read.csv("output/tables/NRI_turnover_ratios.csv", row.names=1)
NTI_turnover_ratios <- read.csv("output/tables/NTI_turnover_ratios.csv", row.names=1)
NRI_turnover_ratios_continental <- read.csv("output/tables/NRI_turnover_ratios_continental.csv", row.names=1)
NTI_turnover_ratios_continental <- read.csv("output/tables/NTI_turnover_ratios_continental.csv", row.names=1)

#download bioclim data
bio <- worldclim_global(var = "bio",res = 10,path = ".")
names(bio)

#temp seasonality (BioClim 4)
temp_seasonality <- bio[[4]]

#annual temperature range (BioClim 7)
annual_temp_range <- bio[[7]]

#extract values for communities
#regional species pool coords
regional_pool_coords <- vect(NRI_turnover_ratios, geom = c("longitude","latitude"), crs = "EPSG:4326")
regional_pool_temp_seasonality_values <- extract(temp_seasonality, regional_pool_coords)
regional_pool_annual_temp_range_values <- extract(annual_temp_range, regional_pool_coords)

#continental species pool coords
continental_pool_coords <- vect(NRI_turnover_ratios_continental, geom = c("longitude","latitude"), crs = "EPSG:4326")
continental_pool_temp_seasonality_values <- extract(temp_seasonality, continental_pool_coords)
continental_pool_annual_temp_range_values <- extract(annual_temp_range, continental_pool_coords)

#add bioclim data to existing dataframes
NRI_turnover_ratios$bio4 <- regional_pool_temp_seasonality_values$wc2.1_10m_bio_4
NTI_turnover_ratios$bio4 <- regional_pool_temp_seasonality_values$wc2.1_10m_bio_4
NRI_turnover_ratios$bio7 <- regional_pool_annual_temp_range_values$wc2.1_10m_bio_7
NTI_turnover_ratios$bio7 <- regional_pool_annual_temp_range_values$wc2.1_10m_bio_7
NRI_turnover_ratios_continental$bio4 <- continental_pool_temp_seasonality_values$wc2.1_10m_bio_4
NTI_turnover_ratios_continental$bio4 <- continental_pool_temp_seasonality_values$wc2.1_10m_bio_4
NRI_turnover_ratios_continental$bio7 <- continental_pool_annual_temp_range_values$wc2.1_10m_bio_7
NTI_turnover_ratios_continental$bio7 <- continental_pool_annual_temp_range_values$wc2.1_10m_bio_7

#plot to check that the values make sense
#example plot
ggplot(NRI_turnover_ratios, aes(abs(latitude), bio4)) +
  geom_point(alpha = 0.3) + geom_smooth(method = "loess") +
  theme_classic() + labs(x = "Latitude", y = "Temperature seasonality (BIO4)")

#plot NRI turnover vs. temperature seasonality (ecoregional) (absolute value of NRI turnover)
NRI_turnover_seasonality_plot <- ggplot(data = NRI_turnover_ratios, aes(x=bio4, y=abs(NRI_turnover))) +
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNRI Difference (June - January) vs. Temperature Seasonality") + 
  labs(x = "Temperature Seasonality (BIO4)", y = "Absolute Value(NRI Difference (June - January))") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE) 

sdmTMB_NRI_turnover_seasonality <- sdmTMB(
  formula = abs(NRI_turnover) ~ bio4,
  data = NRI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NRI_turnover_seasonality, conf.int = TRUE) #while significant, extremely minor slope

#plot NTI turnover vs. temperature seasonality (ecoregional) (NTI turnover)
NTI_turnover_seasonality_plot <- ggplot(data = NTI_turnover_ratios, aes(x=bio4, y=NTI_turnover)) +
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNTI Difference (June - January) vs. Temperature Seasonality") + 
  labs(x = "Temperature Seasonality (BIO4)", y = "NTI Difference (June - January)") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE) 

sdmTMB_NTI_turnover_seasonality <- sdmTMB(
  formula = NTI_turnover ~ bio4,
  data = NTI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NTI_turnover_seasonality, conf.int = TRUE) #more strongly significant negative relationship

### test again with annual temperature range
#plot NRI turnover vs. annual temperature range (ecoregional) (absolute value of NRI turnover)
NRI_turnover_temp_range_plot <- ggplot(data = NRI_turnover_ratios, aes(x=bio7, y=abs(NRI_turnover))) +
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNRI Difference (June - January) vs. Annual Temperature Range") + 
  labs(x = "Annual Temperature Range (BIO7)", y = "Absolute Value (NRI Difference (June - January))") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE) 

sdmTMB_NRI_turnover_temp_range <- sdmTMB(
  formula = abs(NRI_turnover) ~ bio4,
  data = NRI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NRI_turnover_temp_range, conf.int = TRUE) #significant positive slope 
#greater interannual temperature range = greater variation in NRI

#plot NTI turnover vs. temperature temp_range (ecoregional) (NTI turnover)
NTI_turnover_temp_range_plot <- ggplot(data = NTI_turnover_ratios, aes(x=bio4, y=NTI_turnover)) +
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNTI Difference (June - January) vs. Annual Temperature Range") + 
  labs(x = "Annual Temperature Range (BIO7)", y = "NTI Difference (June - January)") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE) 

sdmTMB_NTI_turnover_temp_range <- sdmTMB(
  formula = NTI_turnover ~ bio4,
  data = NTI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NTI_turnover_temp_range, conf.int = TRUE) #more strongly significant negative relationship

#linear models to assess
NRI_seasonality_lm <- lm(abs(NRI_turnover) ~ abs(richness_turnover)*bio4, data=NRI_turnover_ratios)
summary(NRI_seasonality_lm) #greater turnover in species richness and more seasonality predict greater magnitude of NRI turnover 

NTI_seasonality_lm <- lm(abs(NTI_turnover) ~ abs(richness_turnover)*bio4, data=NTI_turnover_ratios)
summary(NTI_seasonality_lm) #greater turnover in species richness and more seasonality predicts lower magnitude of NTI turnover

NRI_temp_range_lm <- lm(abs(NRI_turnover) ~ abs(richness_turnover)*bio7, data=NRI_turnover_ratios)
summary(NRI_temp_range_lm) #greater turnover in species richness and more temp_range predict greater magnitude of NRI turnover 

NTI_temp_range_lm <- lm(abs(NTI_turnover) ~ abs(richness_turnover)*bio7, data=NTI_turnover_ratios)
summary(NTI_temp_range_lm) #greater turnover in species richness and more temp_range predicts lower magnitude of NTI turnover

summary(lm(abs(NRI_turnover) ~ bio4, data=NRI_turnover_ratios)) #positive, p=0.0129
summary(lm(abs(NTI_turnover) ~ bio4, data=NTI_turnover_ratios)) #insignifcant
summary(lm(abs(NRI_turnover) ~ bio7, data=NRI_turnover_ratios)) #positive, p<2e-16
summary(lm(abs(NTI_turnover) ~ bio7, data=NTI_turnover_ratios)) #negative, p=-0.0053

summary(lm(NRI_turnover ~ bio4, data=NRI_turnover_ratios)) #positive, p=2.93e-7
summary(lm(NTI_turnover ~ bio4, data=NTI_turnover_ratios)) #negative, p=2e-16
summary(lm(NRI_turnover ~ bio7, data=NRI_turnover_ratios)) #insignificant
summary(lm(NTI_turnover ~ bio7, data=NTI_turnover_ratios)) #negative, p=2e-16

#this helps us out with NRI/phylogenetic clustering mechanism --
#greater degree of seasonality is allowing for greater magnitude of NRI, which is a larger restructuring 

