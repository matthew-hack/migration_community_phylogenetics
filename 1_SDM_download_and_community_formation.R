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

# define species whose data is to be downloaded from eBirdst
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

#create vector of scientific names of the above species for interfacing with phylogenetic trees downstream
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
