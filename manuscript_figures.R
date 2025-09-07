# 4. Exploratory values and figures
# Generates the results values demonstrating the range of the datasets.
# Generates figures showing range and geographic distributions of NRI and NTI values.
# Generates figures demonstrating variation in species richness and migratory proportion.

#load required libraries
required.libraries <- c("ebirdst","terra","sf","dplyr","rnaturalearth",
                        "ggplot2","viridis","data.table","lme4","lmerTest",
                        "sjPlot")

lapply(required.libraries, require, character.only = TRUE)

project <- terra::project 

# load the mapping parameters to generate necessary CRS
path <- get_species_path("example_data")
fac_parameters <- load_fac_map_parameters(path)
crs <- fac_parameters$custom_projection

#generate gg figure to plot spatial distributions
countries <- ne_countries(returnclass = "sf") %>% 
  st_geometry() %>% 
  st_transform(crs)
country_map <- ggplot()+
  geom_sf(data = countries) +
  coord_sf(xlim=c(-2750000, 3000000), ylim=c(-2000000, 4250000))

#load example raster
ovenbi1_abd_seasonal_mean <- load_raster("data/ebirdst/2021/ovenbi1", product = "abundance", 
                                         period = "seasonal", metric = "mean", 
                                         resolution = "hr")
#use example raster to crop to US boundaries
us <- ne_countries(country = "united states of america", type = "countries", returnclass = "sf") %>% 
  st_transform(st_crs(ovenbi1_abd_seasonal_mean))

##### PLOT SPATIOTEMPORAL PATTERNS OF MIGRATORY PROPORTION (FIGURE 2) ####
#plot migratory proportion distribution for week 1

#load in migratory proportion data for week 1
migratory_proportion_df_week_1 <- read.csv("output/community_data/migratory_proportion_week_1.csv")
migratory_proportion_df_week_1 <- data.frame(migratory_proportion_df_week_1[, c("community", "nmigrant", "nspecies", "proportion")],
                                             row.names = migratory_proportion_df_week_1$community) 

#load first taxonomic species at minimum community size (established in raster_download.R)
wewpew_1_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 1), us), crs, method = "near"), fact=2)

#add coordinates to migratory proportion week 1 data
#generate columns for coordinate data
migratory_proportion_df_week_1.rows <- as.numeric(migratory_proportion_df_week_1$community)
migratory_proportion_df_week_1$x <- rep(NA, nrow(migratory_proportion_df_week_1))
migratory_proportion_df_week_1$y <- rep(NA, nrow(migratory_proportion_df_week_1))
#populate coordinate data
for(i in 1:nrow(migratory_proportion_df_week_1)){
  migratory_proportion_df_week_1$x[i] <- xyFromCell(wewpew_1_agg, migratory_proportion_df_week_1.rows[i])[,1]
  migratory_proportion_df_week_1$y[i] <- xyFromCell(wewpew_1_agg, migratory_proportion_df_week_1.rows[i])[,2]
}

#plot spatial distribution of migratory proportion, week 1 (included in Figure 2)
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
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 30),
        legend.title = element_text(margin = margin(b=20)),
        legend.key.height = unit(2, 'cm')) +
  labs(x = "", y = "")
dev.off()

#plot migratory proportion distribution for week 26
migratory_proportion_df_week_26 <- read.csv("output/community_data/migratory_proportion_week_26.csv")
migratory_proportion_df_week_26 <- data.frame(migratory_proportion_df_week_26[, c("community", "nmigrant", "nspecies", "proportion")],
                                              row.names = migratory_proportion_df_week_26$community) 

#load first taxonomic species at minimum community size (established in raster_download.R)
eawpew_26_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/eawpew", product = "abundance", resolution = "hr"), 26), us), crs, method = "near"), fact=2)

#add coordinates to migratory proportion week 26 data
#generate columns for coordinate data
migratory_proportion_df_week_26.rows <- as.numeric(migratory_proportion_df_week_26$community)
migratory_proportion_df_week_26$x <- rep(NA, nrow(migratory_proportion_df_week_26))
migratory_proportion_df_week_26$y <- rep(NA, nrow(migratory_proportion_df_week_26))
#populate coordinate data
for(i in 1:nrow(migratory_proportion_df_week_26)){
  migratory_proportion_df_week_26$x[i] <- xyFromCell(eawpew_26_agg, migratory_proportion_df_week_26.rows[i])[,1]
  migratory_proportion_df_week_26$y[i] <- xyFromCell(eawpew_26_agg, migratory_proportion_df_week_26.rows[i])[,2]
}

#plot spatial distribution of migratory proportion, week 26 (included in Figure 2)
pdf(file = 'output/figures/community_migratory_proportion_week_26.pdf')
country_map + 
  geom_point(data = migratory_proportion_df_week_26, alpha = 0.1, size = 0.1, aes(x = x, y = y, color = proportion)) +
  scale_color_gradient2("Mig. Prop.", low = "red", high = "blue4",
                        mid = "gray90",
                        midpoint = 0.5) +
  ggtitle("Jun 28") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 30),
        legend.title = element_text(margin = margin(b=20)),
        legend.key.height = unit(2, 'cm')) +
  labs(x = "", y = "")
dev.off()

##### PLOT SPATIOTEMPORAL PATTERNS OF SPECIES RICHNESS (FIGURE 2) ####

#plot spatial distribution of species richness, week 1 (included in Figure 2)
pdf(file = 'output/figures/species_richness_week_1.pdf')
country_map + 
  geom_point(data = migratory_proportion_df_week_1, alpha = 0.1, size = 0.1, aes(x = x, y = y, color = nspecies)) +
  scale_color_viridis("Richness", limits=c(0,92)) +
  ggtitle("Jan 4") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 30),
        legend.title = element_text(margin = margin(b=20)),
        legend.key.height = unit(2, 'cm')) +
  labs(x = "", y = "")
dev.off()

#plot spatial distribution of species richness, week 26 (included in Figure 2)
pdf(file = 'output/figures/species_richness_week_26.pdf')
country_map + 
  geom_point(data = migratory_proportion_df_week_26, alpha = 0.1, size = 0.1, aes(x = x, y = y, color = nspecies)) +
  scale_color_viridis("Richness", limits=c(0,92)) +
  ggtitle("Jun 28") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 30),
        legend.title = element_text(margin = margin(b=20)),
        legend.key.height = unit(2, 'cm')) +
  labs(x = "", y = "")
dev.off()

##### GENERATE HISTOGRAMS OF MIGRATORY PROPORTION AND SPECIES RICHNESS #####

#generate dataframe to populate data for migratory proportion
proportion_hist <- data.frame(matrix(NA, nrow = nrow(migratory_proportion_df_week_26)+nrow(migratory_proportion_df_week_1), ncol = 2))
colnames(proportion_hist) <- c("value", "type")
#populate migratory proportion data, week 1
proportion_hist[c(1:nrow(migratory_proportion_df_week_1)), 1] <- migratory_proportion_df_week_1$proportion
#populate migratory proportion data, week 26
proportion_hist[c((nrow(migratory_proportion_df_week_1) + 1)):(nrow(proportion_hist)),1] <- migratory_proportion_df_week_26$proportion
#add season labels
proportion_hist[c(1:nrow(migratory_proportion_df_week_1)), 2] <- "Non-Breeding"
proportion_hist[c((nrow(migratory_proportion_df_week_1) + 1)):(nrow(proportion_hist)),2] <- "Breeding"

#plot histogram for migratory proportion (summer vs. winter) (included in Figure 2)
pdf(file = 'output/figures/migratory_proportion_histogram.pdf')
ggplot(proportion_hist, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("Proportion Species Migratory") +
  ylab("Density") + 
  guides(fill=guide_legend(title="Season")) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25))
dev.off()

#generate dataframe to populate data for species richness
nspecies_hist <- data.frame(matrix(NA, nrow = nrow(migratory_proportion_df_week_26)+nrow(migratory_proportion_df_week_1), ncol = 2))
colnames(nspecies_hist) <- c("value", "type")
#populate species richness data, week 1
nspecies_hist[c(1:nrow(migratory_proportion_df_week_1)), 1] <- migratory_proportion_df_week_1$nspecies
#populate species richness data, week 26
nspecies_hist[c((nrow(migratory_proportion_df_week_1) + 1)):(nrow(nspecies_hist)),1] <- migratory_proportion_df_week_26$nspecies
#add season labels
nspecies_hist[c(1:nrow(migratory_proportion_df_week_1)), 2] <- "Non-Breeding"
nspecies_hist[c((nrow(migratory_proportion_df_week_1) + 1)):(nrow(nspecies_hist)),2] <- "Breeding"

#plot histogram for species richness (summer vs. winter) (included in Figure 2)
pdf(file = 'output/figures/species_richness_histogram.pdf')
ggplot(nspecies_hist, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("Species Richness") +
  ylab("Density") + 
  guides(fill=guide_legend(title="Season")) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25))
dev.off()

##### GENERATE EXTREME VALUES AND EXPLORATORY CALCULATIONS FOR RESULTS SECTION OF MANUSCRIPT #####

### EXPLORATORY DATA, SPECIES RICHNESS AND MIGRATORY PROPORTION ###

#load complete NRI and NTI datasets generated in 3_community_migratory_proportion.R
all_combined_week_NRIs_df <- read.csv("output/tables/combined_NRI_df.csv")
all_combined_week_NTIs_df <- read.csv("output/tables/combined_NTI_df.csv")
all_combined_week_NRIs_df <- data.table::rbindlist(all_combined_week_NRIs, use.names=TRUE, idcol=TRUE, fill=TRUE) %>% 
  as.data.frame()
all_combined_week_NTIs_df <- data.table::rbindlist(all_combined_week_NTIs, use.names=TRUE, idcol=TRUE, fill=TRUE) %>% 
  as.data.frame()

#calculate exploratory values for species richness
max(na.omit(all_combined_week_NTIs_df$ntaxa)) #109
min(na.omit(all_combined_week_NTIs_df$ntaxa)) #0
median(na.omit(all_combined_week_NTIs_df$ntaxa)) #34
sd(na.omit(all_combined_week_NTIs_df$ntaxa)) #16.13621

#reload migratory proportion data (generated in 3_community_migratory_proportion.R)
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

all_mig_prop_df <- data.table::rbindlist(all_mig_props, use.names=TRUE, idcol=TRUE, fill=TRUE) %>% 
  as.data.frame()

#calculate exploratory values for migratory proportion
max(na.omit(all_combined_week_NTIs_df$ntaxa)) #109
min(na.omit(all_combined_week_NTIs_df$ntaxa)) #0
median(na.omit(all_combined_week_NTIs_df$ntaxa)) #34
sd(na.omit(all_combined_week_NTIs_df$ntaxa)) #16.13621

max(all_mig_prop_df$proportion[is.finite(all_mig_prop_df$proportion)]) #1
min(na.omit(all_mig_prop_df$proportion)) #0
median(na.omit(all_mig_prop_df$proportion)) #0.857
sd(all_mig_prop_df$proportion[is.finite(all_mig_prop_df$proportion)]) #0.124

#calculate exploratory statistics for species richness and migratory proportion (summer)
max(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26,"ntaxa"])) #84
min(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26,"ntaxa"])) #1
mean(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26,"ntaxa"])) #40.08
sd(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26,"ntaxa"])) #11.62
mean(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 26,"proportion"])) #0.885
sd(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 26,"proportion"])) #0.076

#calculate exploratory statistics for species richness and migratory proportion (winter)
max(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1,"ntaxa"])) #78
min(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1,"ntaxa"])) #1
mean(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1,"ntaxa"])) #23.34
sd(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1,"ntaxa"])) #13.85
mean(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 1,"proportion"])) #0.748
sd(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 1,"proportion"])) #0.151

### COMPARING VALUES IN DIFFERENT BCRS ###
#BCR 12
mean(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26 & 
                                         all_combined_week_NTIs_df$BCR == 12,"ntaxa"])) #summer: 58.13
sd(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26 & 
                                       all_combined_week_NTIs_df$BCR == 12,"ntaxa"])) #SD summer: 6.15
mean(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1 & 
                                         all_combined_week_NTIs_df$BCR == 12,"ntaxa"])) #winter: 13.88
sd(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1 & 
                                       all_combined_week_NTIs_df$BCR == 12,"ntaxa"])) #SD winter: 2.47
mean(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1 & 
                                         all_combined_week_NTIs_df$BCR == 12,"ntaxa"])) #winter: 13.88
mean(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 26 & 
                                         all_combined_week_NRIs_df$BCR == 12,"proportion"])) #summer: 0.950
sd(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 26 & 
                                       all_combined_week_NRIs_df$BCR == 12,"proportion"])) #SD summer: 0.017
mean(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 1 & 
                                         all_combined_week_NRIs_df$BCR == 12,"proportion"])) #winter: 0.768
sd(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 1 & 
                                       all_combined_week_NRIs_df$BCR == 12,"proportion"])) #SD winter: 0.079

#BCR 31
mean(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26 & 
                                         all_combined_week_NTIs_df$BCR == 31,"ntaxa"])) #summer: 21.20
sd(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26 & 
                                       all_combined_week_NTIs_df$BCR == 31,"ntaxa"])) #SD summer: 5.24
mean(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1 & 
                                         all_combined_week_NTIs_df$BCR == 31,"ntaxa"])) #winter: 41.10
sd(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1 & 
                                       all_combined_week_NTIs_df$BCR == 31,"ntaxa"])) #SD winter: 6.03
mean(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 26 & 
                                         all_combined_week_NRIs_df$BCR == 31,"proportion"])) #summer: 0.826
sd(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 26 & 
                                       all_combined_week_NRIs_df$BCR == 31,"proportion"])) #SD summer: 0.050
mean(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 1 & 
                                         all_combined_week_NRIs_df$BCR == 31,"proportion"])) #winter: 0.905
sd(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 1 & 
                                       all_combined_week_NRIs_df$BCR == 31,"proportion"])) #SD winter: 0.041

### ECOREGIONAL DATASET EXTREME VALUES ###

#calculate exploratory values for NRI
max(na.omit(all_combined_week_NRIs_df$NRI)) #7.109175
min(na.omit(all_combined_week_NRIs_df$NRI)) #-4.155872
median(na.omit(all_combined_week_NRIs_df$NRI)) #-0.2703786
sd(na.omit(all_combined_week_NRIs_df$NRI)) #1.039971
length(na.omit(all_combined_week_NRIs_df$NRI[all_combined_week_NRIs_df$NRI > 0]))/length(na.omit(all_combined_week_NRIs_df$NRI)) #0.389

#calculate exploratory values for NTI
max(na.omit(all_combined_week_NTIs_df$NTI)) #4.776999
min(na.omit(all_combined_week_NTIs_df$NTI)) #-4.575684
median(na.omit(all_combined_week_NTIs_df$NTI)) #0.3231434
sd(na.omit(all_combined_week_NTIs_df$NTI)) #0.973363
length(na.omit(all_combined_week_NTIs_df$NTI[all_combined_week_NTIs_df$NTI > 0]))/length(na.omit(all_combined_week_NTIs_df$NTI)) #0.626

##### GENERATE VISUALIZATIONS FOR SUPPLEMENTAL FIGURE 2 #####

#generate histograms for NRI and NTI values (ecoregional dataset)

nri_nti_hist <- data.frame(matrix(NA, nrow = nrow(all_combined_week_NRIs_df)+nrow(all_combined_week_NTIs_df), ncol = 2))
colnames(nri_nti_hist) <- c("value", "type")
nri_nti_hist[c(1:nrow(all_combined_week_NRIs_df)), 1] <- all_combined_week_NRIs_df$NRI
nri_nti_hist[c((nrow(all_combined_week_NRIs_df) + 1)):(nrow(nri_nti_hist)),1] <- all_combined_week_NTIs_df$NTI
nri_nti_hist[c(1:nrow(all_combined_week_NRIs_df)), 2] <- "NRI"
nri_nti_hist[c((nrow(all_combined_week_NRIs_df) + 1)):(nrow(nri_nti_hist)),2] <- "NTI"

ggplot(nri_nti_hist, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("NRI vs NTI Distribution (Ecoregion Pool)") +
  ylab("Density") + 
  guides(fill=guide_legend(title="Metric")) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25))

#### convert x and y values to lat and long

all_combined_week_NRIs_df$latitude <- sf_project(from = crs, to = st_crs(4326), all_combined_week_NRIs_df[,c("x","y")])[,2]
all_combined_week_NRIs_df$longitude <- sf_project(from = crs, to = st_crs(4326), all_combined_week_NRIs_df[,c("x","y")])[,1]

all_combined_week_NTIs_df$latitude <- sf_project(from = crs, to = st_crs(4326), all_combined_week_NTIs_df[,c("x","y")])[,2]
all_combined_week_NTIs_df$longitude <- sf_project(from = crs, to = st_crs(4326), all_combined_week_NTIs_df[,c("x","y")])[,1]

#figure out all those positive longitude values
max(all_combined_week_NRIs_df[all_combined_week_NRIs_df$longitude > 0,"y"]) #y = 8672163 -- filter that out


### get extreme values for continental pool (load all values from annual cycle analyses script)

max(na.omit(NRIlistdf$NRI)) #7.109175
min(na.omit(NRIlistdf$NRI)) #-4.155872
median(na.omit(NRIlistdf$NRI)) #-0.2703786
sd(na.omit(NRIlistdf$NRI)) #1.039971
length(na.omit(NRIlistdf$NRI[NRIlistdf$NRI > 0]))/length(na.omit(NRIlistdf$NRI)) #0.389


max(na.omit(NTIlistdf$NTI)) #7.109175
min(na.omit(NTIlistdf$NTI)) #-4.155872
median(na.omit(NTIlistdf$NTI)) #-0.2703786
sd(na.omit(NTIlistdf$NTI)) #1.039971
length(na.omit(NTIlistdf$NTI[NTIlistdf$NTI > 0]))/length(na.omit(NTIlistdf$NTI)) #0.389

#histogram for NRI and NTI values (continental level)
nri_nti_hist_continental <- data.frame(matrix(NA, nrow = nrow(NRIlistdf)+nrow(NTIlistdf), ncol = 2))
colnames(nri_nti_hist_continental) <- c("value", "type")
nri_nti_hist_continental[c(1:nrow(NRIlistdf)), 1] <- NRIlistdf$NRI
nri_nti_hist_continental[c((nrow(NRIlistdf) + 1)):(nrow(nri_nti_hist_continental)),1] <- NTIlistdf$NTI
nri_nti_hist_continental[c(1:nrow(NRIlistdf)), 2] <- "NRI"
nri_nti_hist_continental[c((nrow(NRIlistdf) + 1)):(nrow(nri_nti_hist_continental)),2] <- "NTI"

ggplot(nri_nti_hist_continental, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("NRI vs NTI Distribution (Continental Pool)") +
  ylab("Density") + 
  guides(fill=guide_legend(title="Metric")) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25))
