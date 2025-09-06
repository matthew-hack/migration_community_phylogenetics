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

# load list of species names to be included in dataset
species_names <- c("rotbec", "nobtyr", "olsfly", "grepew", "wewpew", "eawpew", "yebfly", "acafly",
                   "aldfly", "wilfly", "leafly", "hamfly", "gryfly", "dusfly", "pasfly", "corfly", 
                   "bubfly", "blkpho", "easpho", "saypho", "verfly", "ducfly", "astfly", "grcfly", 
                   "bncfly", "grekis", "subfly", "trokin", "coukin", "caskin", "thbkin", "weskin", 
                   "easkin", "grykin", "sctfly", "bkcvir1", "whevir", "belvir", "gryvir", "hutvir", 
                   "yetvir", "casvir", "buhvir", "plsvir", "phivir", "warvir", "reevir1", "bkwvir", 
                   "logshr", "gryjay", "brnjay", "grnjay", "pinjay", "stejay", "blujay", "flsjay", 
                   "issjay", "cowscj1", "bkbmag1", "yebmag", "clanut", "amecro", "fiscro", "chirav", 
                   "comrav", "carchi", "bkcchi", "mouchi", "mexchi", "chbchi", "borchi2", "gyhchi", 
                   "britit", "oaktit", "juntit1", "tuftit", "blctit4", "verdin", "horlar", "nrwswa", 
                   "purmar", "treswa", "vigswa", "banswa", "cliswa", "cavswa", "rewbul", "arcwar1",
                   "bushti", "wrenti", "ruckin", "gockin", "rebnut", "whbnut", "pygnut", "bnhnut", 
                   "brncre", "buggna", "bktgna", "calgna", "bkcgna", "rocwre", "canwre", "houwre", 
                   "winwre3", "marwre", "carwre", "bewwre", "cacwre", "amedip", "eursta", "commyn", 
                   "grycat", "cubthr", "brnthr", "lobthr", "benthr", "calthr", "lecthr", "crithr", 
                   "sagthr", "normoc", "easblu", "wesblu", "moublu", "towsol", "veery", "gycthr", 
                   "bicthr", "swathr", "herthr", "woothr", "clcrob", "amerob", "norwhe", "bohwax", 
                   "cedwax", "phaino", "oliwar", "orabis1", "nutman", "houspa", "eutspa", "eaywag", 
                   "amepip", "sprpip", "evegro", "pingro", "gcrfin", "bkrfin", "bcrfin", "houfin", 
                   "purfin", "casfin","comred", "hoared", "redcro", "whwcro", "pinsis", "lesgol", 
                   "lawgol", "amegfi", "laplon", "chclon", "smilon", "mcclon", "snobun", "mckbun", 
                   "ruwspa", "botspa", "casspa", "bacspa", "graspa", "olispa", "chispa", "clcspa", 
                   "bkcspa", "fiespa", "brespa", "fisspa", "bktspa", "larspa", "larbun", "amtspa", 
                   "foxspa", "daejun", "yeejun", "whcspa", "gocspa", "harspa", "whtspa", "belspa2", 
                   "vesspa", "lecspa", "seaspa", "nstspa", "sstspa", "savspa", "baispa", "henspa", 
                   "sonspa", "linspa", "swaspa", "cantow", "abetow", "caltow", "rucspa", "gnttow", 
                   "spotow", "eastow", "yebcha", "yehbla", "boboli", "wesmea", "easmea", "orcori", 
                   "hooori", "bulori", "spbori", "altori", "audori", "balori", "scoori", "rewbla", 
                   "tribla", "shicow", "brocow", "bnhcow", "rusbla", "brebla", "comgra", "botgra", 
                   "grtgra", "ovenbi1", "woewar1", "louwat", "norwat", "gowwar", "bawwar", "prowar",
                   "swawar", "tenwar", "orcwar", "colwar", "lucwar", "naswar", "virwar", "conwar", 
                   "macwar", "mouwar", "kenwar", "comyel", "hoowar", "amered", "kirwar", "camwar", 
                   "cerwar", "norpar", "tropar", "magwar", "babwar", "bkbwar", "yelwar", "chswar", 
                   "bkpwar", "btbwar", "palwar", "pinwar", "yerwar", "yetwar", "prawar", "grawar", 
                   "btywar", "towwar", "herwar", "gchwar", "btnwar", "canwar", "wlswar", "refwar", 
                   "paired", "heptan", "sumtan", "scatan", "westan", "norcar", "pyrrhu", "robgro", 
                   "bkhgro", "blugrb1", "lazbun", "indbun", "varbun", "paibun", "dickci")

#create vector of scientific names of the above species for interfacing with phylogenetic tree
species_passerines <- c("Pachyramphus_aglaiae", "Camptostoma_imberbe","Contopus_cooperi",
                        "Contopus_pertinax", "Contopus_sordidulus", "Contopus_virens",
                        "Empidonax_flaviventris", "Empidonax_virescens", "Empidonax_alnorum",
                        "Empidonax_traillii", "Empidonax_minimus", "Empidonax_hammondii",
                        "Empidonax_wrightii", "Empidonax_oberholseri", "Empidonax_difficilis",
                        "Empidonax_occidentalis", "Empidonax_fulvifrons", "Sayornis_nigricans",
                        "Sayornis_phoebe", "Sayornis_saya", "Pyrocephalus_rubinus",
                        "Myiarchus_tuberculifer", "Myiarchus_cinerascens", "Myiarchus_crinitus",
                        "Myiarchus_tyrannulus", "Pitangus_sulphuratus", "Myiodynastes_luteiventris",
                        "Tyrannus_melancholicus", "Tyrannus_couchii", "Tyrannus_vociferans",
                        "Tyrannus_crassirostris", "Tyrannus_verticalis", "Tyrannus_tyrannus",
                        "Tyrannus_dominicensis", "Tyrannus_forficatus", "Vireo_atricapilla",
                        "Vireo_griseus", "Vireo_bellii", "Vireo_vicinior",
                        "Vireo_huttoni", "Vireo_flavifrons", "Vireo_cassinii",
                        "Vireo_solitarius", "Vireo_plumbeus", "Vireo_philadelphicus",
                        "Vireo_gilvus", "Vireo_olivaceus", "Vireo_altiloquus",
                        "Lanius_ludovicianus", "Perisoreus_canadensis", "Cyanocorax_morio",
                        "Cyanocorax_yncas", "Gymnorhinus_cyanocephalus", "Cyanocitta_stelleri",
                        "Cyanocitta_cristata", "Aphelocoma_coerulescens", "Aphelocoma_insularis",
                        "Aphelocoma_californica", "Pica_hudsonia", "Pica_nuttalli",
                        "Nucifraga_columbiana", "Corvus_brachyrhynchos", "Corvus_ossifragus", 
                        "Corvus_cryptoleucus", "Corvus_corax", "Parus_carolinensis", 
                        "Parus_atricapillus", "Parus_gambeli", "Parus_sclateri", 
                        "Parus_rufescens", "Parus_hudsonicus", "Parus_cinctus", 
                        "Baeolophus_wollweberi", "Baeolophus_inornatus", "Baeolophus_ridgwayi",
                        "Baeolophus_bicolor", "Baeolophus_atricristatus", "Auriparus_flaviceps", 
                        "Eremophila_alpestris", "Stelgidopteryx_serripennis", "Progne_subis",
                        "Tachycineta_bicolor", "Tachycineta_thalassina", "Riparia_riparia",
                        "Petrochelidon_pyrrhonota", "Petrochelidon_fulva", "Pycnonotus_jocosus", 
                        "Phylloscopus_borealis", "Psaltriparus_minimus", "Chamaea_fasciata",
                        "Regulus_calendula", "Regulus_satrapa", "Sitta_canadensis",
                        "Sitta_carolinensis", "Sitta_pygmaea", "Sitta_pusilla",
                        "Certhia_americana", "Polioptila_caerulea", "Polioptila_melanura",
                        "Polioptila_californica", "Polioptila_nigriceps", "Salpinctes_obsoletus",
                        "Catherpes_mexicanus", "Troglodytes_aedon", "Troglodytes_troglodytes",
                        "Cistothorus_palustris", "Thryothorus_ludovicianus", "Thryomanes_bewickii",
                        "Campylorhynchus_brunneicapillus", "Cinclus_mexicanus", "Sturnus_vulgaris",
                        "Acridotheres_tristis", "Dumetella_carolinensis", "Toxostoma_curvirostre",
                        "Toxostoma_rufum", "Toxostoma_longirostre", "Toxostoma_bendirei",
                        "Toxostoma_redivivum", "Toxostoma_lecontei", "Toxostoma_crissale",
                        "Oreoscoptes_montanus", "Mimus_polyglottos", "Sialia_sialis",
                        "Sialia_mexicana", "Sialia_currucoides", "Myadestes_townsendi",
                        "Catharus_fuscescens", "Catharus_minimus", "Catharus_bicknelli",
                        "Catharus_ustulatus", "Catharus_guttatus", "Hylocichla_mustelina",
                        "Turdus_grayi", "Turdus_migratorius",
                        "Oenanthe_oenanthe", "Bombycilla_garrulus", "Bombycilla_cedrorum",
                        "Phainopepla_nitens", "Peucedramus_taeniatus","Euplectes_franciscanus",
                        "Lonchura_punctulata", "Passer_domesticus", "Passer_montanus",
                        "Motacilla_flava", "Anthus_rubescens", "Anthus_spragueii",
                        "Coccothraustes_vespertinus", "Pinicola_enucleator", "Leucosticte_tephrocotis",
                        "Leucosticte_atrata", "Leucosticte_australis", "Carpodacus_mexicanus",
                        "Carpodacus_purpureus", "Carpodacus_cassinii", "Carduelis_flammea",
                        "Carduelis_hornemanni", "Loxia_curvirostra", "Loxia_leucoptera",
                        "Carduelis_pinus", "Carduelis_psaltria", "Carduelis_lawrencei",
                        "Carduelis_tristis", "Calcarius_lapponicus", "Calcarius_ornatus",
                        "Calcarius_pictus", "Calcarius_mccownii", "Plectrophenax_nivalis",
                        "Plectrophenax_hyperboreus","Aimophila_carpalis", "Aimophila_botterii",
                        "Aimophila_cassinii", "Aimophila_aestivalis", "Ammodramus_savannarum",
                        "Arremonops_rufivirgatus", "Spizella_passerina", "Spizella_pallida",
                        "Spizella_atrogularis", "Spizella_pusilla", "Spizella_breweri",
                        "Aimophila_quinquestriata", "Amphispiza_bilineata", "Chondestes_grammacus",
                        "Calamospiza_melanocorys", "Spizella_arborea", "Passerella_iliaca",
                        "Junco_hyemalis", "Junco_phaeonotus", "Zonotrichia_leucophrys",
                        "Zonotrichia_atricapilla", "Zonotrichia_querula", "Zonotrichia_albicollis",
                        "Amphispiza_belli", "Pooecetes_gramineus", "Ammodramus_leconteii",
                        "Ammodramus_maritimus", "Ammodramus_nelsoni", "Ammodramus_caudacutus",
                        "Passerculus_sandwichensis", "Ammodramus_bairdii", "Ammodramus_henslowii",
                        "Melospiza_melodia", "Melospiza_lincolnii", "Melospiza_georgiana",
                        "Pipilo_fuscus", "Pipilo_aberti", "Pipilo_crissalis",
                        "Aimophila_ruficeps", "Pipilo_chlorurus", "Pipilo_maculatus",
                        "Pipilo_erythrophthalmus", "Icteria_virens", "Xanthocephalus_xanthocephalus",
                        "Dolichonyx_oryzivorus", "Sturnella_neglecta", "Sturnella_magna",
                        "Icterus_spurius", "Icterus_cucullatus", "Icterus_bullockii",
                        "Icterus_pectoralis", "Icterus_gularis", "Icterus_graduacauda",
                        "Icterus_galbula", "Icterus_parisorum", "Agelaius_phoeniceus",
                        "Agelaius_tricolor", "Molothrus_bonariensis", "Molothrus_aeneus",
                        "Molothrus_ater", "Euphagus_carolinus", "Euphagus_cyanocephalus",
                        "Quiscalus_quiscula", "Quiscalus_major", "Quiscalus_mexicanus",
                        "Seiurus_aurocapilla", "Helmitheros_vermivorum", "Seiurus_motacilla",
                        "Seiurus_noveboracensis", "Vermivora_chrysoptera", "Mniotilta_varia",
                        "Protonotaria_citrea", "Limnothlypis_swainsonii", "Vermivora_peregrina",
                        "Vermivora_celata", "Vermivora_crissalis", "Vermivora_luciae",
                        "Vermivora_ruficapilla", "Vermivora_virginiae", "Oporornis_agilis",
                        "Oporornis_tolmiei", "Oporornis_philadelphia", "Oporornis_formosus",
                        "Geothlypis_trichas", "Wilsonia_citrina", "Setophaga_ruticilla",
                        "Dendroica_kirtlandii", "Dendroica_tigrina", "Dendroica_cerulea",
                        "Parula_americana", "Parula_pitiayumi", "Dendroica_magnolia",
                        "Dendroica_castanea", "Dendroica_fusca", "Dendroica_petechia",
                        "Dendroica_pensylvanica", "Dendroica_striata", "Dendroica_caerulescens",
                        "Dendroica_palmarum", "Dendroica_pinus", "Dendroica_coronata",
                        "Dendroica_dominica", "Dendroica_discolor", "Dendroica_graciae",
                        "Dendroica_nigrescens", "Dendroica_townsendi", "Dendroica_occidentalis",
                        "Dendroica_chrysoparia", "Dendroica_virens", "Wilsonia_canadensis",
                        "Wilsonia_pusilla", "Cardellina_rubrifrons", "Myioborus_pictus",
                        "Piranga_flava", "Piranga_rubra", "Piranga_olivacea",
                        "Piranga_ludoviciana", "Cardinalis_cardinalis", "Cardinalis_sinuatus",
                        "Pheucticus_ludovicianus", "Pheucticus_melanocephalus", "Passerina_caerulea",
                        "Passerina_amoena", "Passerina_cyanea", "Passerina_versicolor",
                        "Passerina_ciris", "Spiza_americana")

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

#apply BCR values to NRI list (warning! takes long processing time)
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

#export dataframe
write.csv(NRIlistdf, "output/tables/all_NRI_data.csv")

#reload dataframe if already processed
NRIlistdf <- read.csv("output/tables/all_NRI_data.csv")


##### ADD MIGRATORY PROPORTION TO THE OTHER WEEKS #####

for(i in 1:52){ #loading in the files takes about a half hour
  assign(paste0("migratory_proportion_df_week", i), read.csv(paste0("output/community_data/migratory_proportion_week_", i, ".csv")))
} 

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

all_combined_week_NRIs_df$proportion <- rep(NA, nrow(all_combined_week_NRIs_df))
for(i in 1:nrow(all_combined_week_NRIs_df)){ 
  all_combined_week_NRIs_df[i,"proportion"] <- all_mig_props[[all_combined_week_NRIs_df[i,".id"]]][all_mig_props[[all_combined_week_NRIs_df[i,".id"]]]$community == all_combined_week_NRIs_df[i,"X"],"proportion"]
}
#warning: takes a long time (~3 days)
write.csv(all_combined_week_NRIs_df, "output/tables/combined_NRI_df.csv")
all_combined_week_NRIs_df <- read.csv("output/tables/combined_NRI_df.csv")

#replicate with NTI
all_combined_week_NTIs_df$proportion <- rep(NA, nrow(all_combined_week_NTIs_df))
for(i in 1:nrow(all_combined_week_NTIs_df)){ 
  all_combined_week_NTIs_df[i,"proportion"] <- all_mig_props[[all_combined_week_NTIs_df[i,".id"]]][all_mig_props[[all_combined_week_NTIs_df[i,".id"]]]$community == all_combined_week_NTIs_df[i,"X"],"proportion"]
}
#warning: takes a long time (~3 days)
write.csv(all_combined_week_NTIs_df, "output/tables/combined_NTI_df.csv")
all_combined_week_NTIs_df <- read.csv("output/tables/combined_NTI_df.csv")
