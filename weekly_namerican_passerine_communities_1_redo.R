# 2. Weekly community passerine phylogenetic metric generation
# This script runs all community phylogenetic metrics for one week of the year (week 1)
# NRI and NTI metrics generated for each community at two scales of species pools
# Creates visualizations for spatial distribution of NRI and NTI metrics at both species pools
# Script can be replicated for other 51 weeks of the year with minimal substitutions (indicated below) 

#load required libraries
required.libraries <- c("ebirdst","terra","sf","dplyr","tidyr","rnaturalearth",
                        "geodata","ggplot2","fields","stringr","TreeTools",
                        "viridis","rgdal","raster","devtools","data.table",
                        "ape","picante","pez","phytools","vegan","adephylo",
                        "phylobase","geiger","mvMORPH","OUwie","hisse",
                        "BAMMtools","phylosignal","Biostrings","kableExtra",
                        "betapart","gridExtra","reshape2","data.table","ggmap")

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

##### PREPARE COMMUNITY AND PHYLOGENY FILES FOR COMMUNITY PHYLOGENETICS METRICS #####

#load first taxonomic species at minimum community size (established in raster_download.R)
wewpew_1_agg_df <- read.csv("output/species_data/1/wewpew_1_agg.csv", row.names = 1)
#read existing community table (generated in raster_download.R)
community1 <- data.table::fread("output/community_data/community_1.csv")
#assign rownames to match the individually identified communities
rownames(community1) <- rownames(wewpew_1_agg_df)
#remove "V1" columns if they have been generated in loading in community data
community1 <- data.frame(community1[,-c("V1")])

#clean phylogenetic trees and community matrix using picante functions
#clean phylogenetic tree
passerine.picCleanTree <- match.phylo.comm(phy = passerine.tree, 
                                           comm = community1)$phy
#export tree
write.tree(passerine.picCleanTree, file = "passerine.picCleanTree")
#clean community matrix
passerine.picCleanComm <- match.phylo.comm(phy = passerine.tree, 
                                           comm = community1)$comm

##### RUN COMMUNITY PHYLOGENETIC METRICS ON CONTINENTAL SPECIES POOL #####
#names and file structure specific to community 1 -- each community data has its own folder in the "output/tables" file path

#subsample dataframe
set.seed(52)
#create vector of positions in community file
ints <- seq(1, nrow(community1), 1) 
#subsample 25000 communities for computatinal purposes
intsample <- sample(ints, 25000)
subsampled_communities <- as.numeric(row.names(passerine.picCleanComm[intsample,]))
write.csv(as.data.frame(subsampled_communities), "output/tables/1/subsampled_communities.csv")
#reload subsampled communities if already created
intsample_df <- read.csv("output/tables/1/subsampled_communities.csv", row.names=1)
intsample <- unlist(as.vector(intsample_df))

#randomization -- SES PD null model
passerine.1.picCleanComm.subsample <- passerine.picCleanComm[rownames(passerine.picCleanComm) %in% as.character(intsample), ] 

#create a phylogenetic distance matrix
passerine.cophenDist <- cophenetic.phylo(passerine.picCleanTree)

#generate MPD (mean phylogenetic distance) values for each community at continental species pool
passerine.ses.mpd.1 <- ses.mpd(passerine.1.picCleanComm.subsample,
                                passerine.cophenDist,
                                null.model = "taxa.labels",
                                abundance.weighted = FALSE,
                                runs = 100)
#export and reread data frame containing MPD values
write.csv(passerine.ses.mpd.1, "output/tables/1/passerine.ses.mpd.1.csv")
passerine.ses.mpd.1 <- read.csv("output/tables/1/passerine.ses.mpd.1.csv", row.names=1, header=TRUE)

#generate MNTD (mean nearest taxa distance) values for each community at continental species pool
passerine.ses.mntd.1 <- ses.mntd(passerine.1.picCleanComm.subsample,
                                  passerine.cophenDist,
                                  null.model = "taxa.labels",
                                  abundance.weighted = FALSE,
                                  runs = 100)
#export and reread data frame containing MNTD values
write.csv(passerine.ses.mntd.1, "output/tables/1/passerine.ses.mntd.1.csv")
passerine.ses.mntd.1 <- read.csv("output/tables/1/passerine.ses.mntd.1.csv", row.names=1, header=TRUE)

# visualizations

#passerine.1.ses.pd.subsample.rows <- as.numeric(row.names(passerine.1.ses.pd.subsample))

#passerine.1.ses.pd.subsample$x <- rep(NA, nrow(passerine.1.ses.pd.subsample))
#passerine.1.ses.pd.subsample$y <- rep(NA, nrow(passerine.1.ses.pd.subsample))

#use agg spatraster of minimum length df
#wewpew_1_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 1), us), crs, method = "near"), fact=2)

#add coords to each subsampled community
#for(i in 1:nrow(passerine.1.ses.pd.subsample)){
#  passerine.1.ses.pd.subsample$x[i] <- xyFromCell(wewpew_1_agg, passerine.1.ses.pd.subsample.rows[i])[,1]
#  passerine.1.ses.pd.subsample$y[i] <- xyFromCell(wewpew_1_agg, passerine.1.ses.pd.subsample.rows[i])[,2]
#}

#plot community SES PD
#ggplot(data = passerine.1.ses.pd.subsample, mapping = aes(x = x, y = y)) +
#  geom_point(alpha = 0.5, aes(color = pd.obs.z)) +
#  scale_color_gradient2("SES PD", low = "red", high = "blue4",
#                        mid = "gray90",
#                        midpoint = 0, limits = c(min(passerine.1.ses.pd.subsample$pd.obs.z), max(passerine.1.ses.pd.subsample$pd.obs.z)),
#                       na.value = "gray97") +
#  theme(axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank(),
#        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        plot.title = element_text(hjust = 0.5)) +
#  labs(x = "", y = "")

#plot community NRI
#ggplot(data = passerine.1.ses.pd.subsample, mapping = aes(x = x, y = y)) + 
#  geom_point(alpha = 0.5, aes(color = pd.obs.z)) + 
#  scale_color_gradient2("SES PD", low = "red", high = "blue4", 
#                        mid = "gray90",
#                        midpoint = 0, limits = c(min(passerine.1.ses.pd.subsample$pd.obs.z), max(passerine.1.ses.pd.subsample$pd.obs.z)),
#                        na.value = "gray97") +  
#  theme(axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank(),
#        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        plot.title = element_text(hjust = 0.5)) +
#  labs(x = "", y = "")



#Visualizations

# 
# passerine.1.ses.pd.subsample <- read.csv("output/tables/1/passerine.pd.1.subsample.csv", row.names=1, header=TRUE)
# passerine.1.ses.pd.subsample.rows <- as.numeric(row.names(passerine.1.ses.pd.subsample))
# 
# passerine.1.ses.pd.subsample$x <- rep(NA, nrow(passerine.1.ses.pd.subsample))
# passerine.1.ses.pd.subsample$y <- rep(NA, nrow(passerine.1.ses.pd.subsample))
# 
# #use agg spatraster of minimum length df
# wewpew_1_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 1), us), crs, method = "near"), fact=2)
# 
# #add coords to each subsampled community
# for(i in 1:nrow(passerine.1.ses.pd.subsample)){
#   passerine.1.ses.pd.subsample$x[i] <- xyFromCell(wewpew_1_agg, passerine.1.ses.pd.subsample.rows[i])[,1]
#   passerine.1.ses.pd.subsample$y[i] <- xyFromCell(wewpew_1_agg, passerine.1.ses.pd.subsample.rows[i])[,2]
# }
# 
# #plot community SES PD
# #pdf(file = 'output/figures/1/passerine.1.PD.pdf')
# ggplot(data = passerine.1.ses.pd.subsample, mapping = aes(x = x, y = y)) +
#   geom_point(alpha = 0.5, aes(color = pd.obs.z)) +
#   scale_color_gradient2("SES PD", low = "blue4", high = "red",
#                         mid = "gray90",
#                         midpoint = 0, limits = c(min(passerine.1.ses.pd.subsample$pd.obs.z), max(passerine.1.ses.pd.subsample$pd.obs.z)),
#                         na.value = "gray97") +
#   ggtitle("Jan 4") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5)) +
#   labs(x = "", y = "")
# #dev.off()
# # 
# # 
# # 
# # #NRI
# passerine.ses.mpd.1 <- read.csv("output/tables/1/passerine.ses.mpd.1.csv", row.names=1, header=TRUE)
# passerine.ses.mpd.1.rows <- as.numeric(row.names(passerine.ses.mpd.1))
# 
# passerine.ses.mpd.1$x <- rep(NA, nrow(passerine.ses.mpd.1))
# passerine.ses.mpd.1$y <- rep(NA, nrow(passerine.ses.mpd.1))
# passerine.ses.mpd.1$NRI <- rep(NA, nrow(passerine.ses.mpd.1))
# 
# 
# #add coords to each subsampled community
# for(i in 1:nrow(passerine.ses.mpd.1)){
#   passerine.ses.mpd.1$x[i] <- xyFromCell(wewpew_1_agg, passerine.ses.mpd.1.rows[i])[,1]
#   passerine.ses.mpd.1$y[i] <- xyFromCell(wewpew_1_agg, passerine.ses.mpd.1.rows[i])[,2]
# }
# 
# passerine.ses.mpd.1$NRI <- -1*passerine.ses.mpd.1$mpd.obs.z
# 
# library(ggmap)
# 
# #plot community NRI
# country_map <- ggplot()+
#   geom_sf(data = countries) +
#   coord_sf(xlim=c(-2750000, 3000000), ylim=c(-2000000, 4250000))
# 
# pdf(file = 'output/figures/1/passerine.1.NRI.pdf')
# country_map + 
#   geom_point(data = passerine.ses.mpd.1, alpha = 0.5, aes(x = x, y = y, color = NRI)) +
#   scale_color_gradient2("NRI", low = "red", high = "blue4",
#                         mid = "gray90",
#                         midpoint = 0, limits = c(-4,7),
#                         na.value = "gray97") +
#   ggtitle("Jan 4") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5)) +
#   labs(x = "", y = "")
# dev.off()
# 
# country_map + 
#   geom_point(data = passerine.ses.mpd.1, alpha = 0.5, aes(x = x, y = y, color = ntaxa)) +
#   scale_color_gradient2("Richness", low = "red", high = "blue4",
#                         mid = "gray90",
#                         midpoint = median(passerine.ses.mpd.1$ntaxa)) +
#   ggtitle("Jan 4") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5)) +
#   labs(x = "", y = "")
# 
# head(passerine.ses.mpd.1)
# 
# 
# #NTI
# passerine.ses.mntd.1 <- read.csv("output/tables/1/passerine.ses.mntd.1.csv", row.names=1, header=TRUE)
# passerine.ses.mntd.1.rows <- as.numeric(row.names(passerine.ses.mntd.1))
# 
# passerine.ses.mntd.1$x <- rep(NA, nrow(passerine.ses.mntd.1))
# passerine.ses.mntd.1$y <- rep(NA, nrow(passerine.ses.mntd.1))
# passerine.ses.mntd.1$NTI <- rep(NA, nrow(passerine.ses.mntd.1))
# 
# 
# #add coords to each subsampled community
# for(i in 1:nrow(passerine.ses.mntd.1)){
#   passerine.ses.mntd.1$x[i] <- xyFromCell(wewpew_1_agg, passerine.ses.mntd.1.rows[i])[,1]
#   passerine.ses.mntd.1$y[i] <- xyFromCell(wewpew_1_agg, passerine.ses.mntd.1.rows[i])[,2]
# }
# 
# passerine.ses.mntd.1$NTI <- -1*passerine.ses.mntd.1$mntd.obs.z
# 
# 
# #plot community NTI
# country_map <- ggplot()+
#   geom_sf(data = countries) +
#   coord_sf(xlim=c(-2750000, 3000000), ylim=c(-2000000, 4250000))
# 
# pdf(file = 'output/figures/1/passerine.1.NTI.pdf')
# country_map + 
#   geom_point(data = passerine.ses.mntd.1, alpha = 0.5, aes(x = x, y = y, color = NTI)) +
#   scale_color_gradient2("NTI", low = "red", high = "blue4",
#                         mid = "gray90",
#                         midpoint = 0, limits = c(-4,7),
#                         na.value = "gray97") +
#   ggtitle("Jan 4") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5)) +
#   labs(x = "", y = "")
# dev.off()
# 
# 
# ### Attempt BCR-level communities
# 
# # bcr12 <- read_sf("data/gis-data.gpkg", "bcr12") %>% 
# #   st_transform(crs = crs) %>% 
# #   st_geometry()
# # wewpew_1_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 1), us), crs, method = "near"), fact=2)
# # wewpew_1_bcr_12 <- rast(mask(raster(wewpew_1_agg), as_Spatial(bcr12)))
# # wewpew_1_bcr_12_df <- as.data.frame(wewpew_1_bcr_12)
# # bcr12rows <- row.names(wewpew_1_bcr_12_df)
# # 
# # community_1_bcr_12 <- community1[bcr12rows,]
# # community_1_bcr_12_fil <- community_1_bcr_12[, colSums(community_1_bcr_12 != 0) > 0] #looks wrong btw
# # community_1_bcr_12_fil
# # 
# # passerine.picCleanTree_bcr12 <- match.phylo.comm(phy = passerine.tree, 
# #                                                  comm = community_1_bcr_12_fil)$phy
# # passerine.picCleanComm_bcr12 <- match.phylo.comm(phy = passerine.tree, 
# #                                                  comm = community_1_bcr_12_fil)$comm
# # 
# # ints12 <- seq(1, nrow(community_1_bcr_12_fil), 1) #create vector of positions
# # intsample12 <- sample(ints12, 10000) #sample positions
# # 
# # #randomization -- SES PD null model
# # passerine.1.picCleanComm.subsample.bcr12 <- passerine.picCleanComm_bcr12[intsample12,] #feel free to modify
# # passerine.cophenDist.bcr12 <- cophenetic.phylo(passerine.picCleanTree_bcr12)
# # 
# # #NRI on summer community
# # passerine.ses.mpd.1.bcr12 <- ses.mpd(passerine.1.picCleanComm.subsample.bcr12,
# #                                       passerine.cophenDist.bcr12,
# #                                       null.model = "taxa.labels",
# #                                       abundance.weighted = FALSE,
# #                                       runs = 100)
# # #NTI on summer community
# # passerine.ses.mntd.1.bcr12 <- ses.mntd(passerine.1.picCleanComm.subsample.bcr12,
# #                                         passerine.cophenDist.bcr12,
# #                                         null.model = "taxa.labels",
# #                                         abundance.weighted = FALSE,
# #                                         runs = 100)
# # 
# # 
# # 
# # #visualize
# # passerine.ses.mpd.1.bcr12.rows <- as.numeric(row.names(passerine.ses.mpd.1.bcr12))
# # 
# # passerine.ses.mpd.1.bcr12$x <- rep(NA, nrow(passerine.ses.mpd.1.bcr12))
# # passerine.ses.mpd.1.bcr12$y <- rep(NA, nrow(passerine.ses.mpd.1.bcr12))
# # passerine.ses.mpd.1.bcr12$NRI <- rep(NA, nrow(passerine.ses.mpd.1.bcr12))
# # 
# # 
# # #add coords to each subsampled community
# # for(i in 1:nrow(passerine.ses.mpd.1.bcr12)){
# #   passerine.ses.mpd.1.bcr12$x[i] <- xyFromCell(wewpew_1_bcr_12, passerine.ses.mpd.1.bcr12.rows[i])[,1]
# #   passerine.ses.mpd.1.bcr12$y[i] <- xyFromCell(wewpew_1_bcr_12, passerine.ses.mpd.1.bcr12.rows[i])[,2]
# # }
# # 
# # passerine.ses.mpd.1.bcr12$NRI <- -1*passerine.ses.mpd.1.bcr12$mpd.obs.z
# # 
# # 
# # #plot community NRI
# # ggplot(data = passerine.ses.mpd.1.bcr12, mapping = aes(x = x, y = y)) +
# #   geom_point(alpha = 0.5, aes(color = NRI)) +
# #   scale_color_gradient2("NRI", low = "red", high = "blue4",
# #                         mid = "gray90",
# #                         midpoint = 0, limits = c(min(passerine.ses.mpd.1.bcr12$NRI), max(passerine.ses.mpd.1.bcr12$NRI)),
# #                         na.value = "gray97") +
# #   ggtitle("Jan 4") +
# #   theme(axis.text.x=element_blank(),
# #         axis.ticks.x=element_blank(),
# #         axis.text.y=element_blank(),
# #         axis.ticks.y=element_blank(),
# #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
# #         plot.title = element_text(hjust = 0.5)) +
# #   labs(x = "", y = "")
# # 
# # ggplot(data = passerine.ses.mpd.1.bcr12, mapping = aes(x = x, y = y)) +
# #   geom_point(alpha = 0.5, aes(color = ntaxa)) +
# #   scale_color_viridis("Richness") +
# #   ggtitle("Jan 4") +
# #   theme(axis.text.x=element_blank(),
# #         axis.ticks.x=element_blank(),
# #         axis.text.y=element_blank(),
# #         axis.ticks.y=element_blank(),
# #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
# #         plot.title = element_text(hjust = 0.5)) +
# #   labs(x = "", y = "")
# # 
# # ggplot(data = passerine.ses.mpd.1.bcr12, mapping = aes(x=ntaxa, y=NRI)) +
# #   geom_point(alpha = 0.1)
# 
# # #load in minimum length species vector
wewpew_1_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 1), us), crs, method = "near"), fact=2)
# #read in all BCR data
allbcrs <- read_sf("data/BCR_Terrestrial/BCR_Terrestrial_master.shp") #read in BCR files
#assign row names to BCRs
for(i in 1:37){
  assign(paste0("bcr",i), allbcrs[allbcrs$BCR == i,] %>%
           st_transform(crs = crs) %>%
           st_geometry())
  assign(paste0("bcr", i, "spatial"), as_Spatial(get(paste0("bcr", i))))
  assign(paste0("wewpew_1_bcr_", i), rast(mask(raster(wewpew_1_agg), get(paste0("bcr", i, "spatial")))))
  assign(paste0("wewpew_1_bcr_", i, "_df"), as.data.frame(get(paste0("wewpew_1_bcr_", i))))
  assign(paste0("bcr", i, "rows"), row.names(get(paste0("wewpew_1_bcr_", i, "_df"))))
}

# #read in community data
community1 <- data.table::fread("output/community_data/community_1.csv")
#deals with the issue of loading in multiple V1 copies but also preserves row names
community1 <- data.frame(community1[,!"V1"], row.names = community1$V1)

#create weekly community data subsetted by BCR
for(i in 1:37){
  assign(paste0("community_1_bcr_", i), subset(community1, row.names(community1) %in% get(paste0("bcr", i, "rows"))))
  assign(paste0("community_1_bcr_", i, "_fil"),
         get(paste0("community_1_bcr_", i)) %>% dplyr::select(which(colSums(get(paste0("community_1_bcr_", i))) > 0)))
}
#export the data
for(i in 1:37){
  data.table::fwrite(get(paste0("community_1_bcr_", i, "_fil")), file = paste0("output/bcr_community_data/1/community_1_bcr_", i, ".csv"))
}

# Using picante -- clean trees/comm
for(i in 1:37){
  assign(paste0("passerine.picCleanTree.", i), match.phylo.comm(phy = passerine.tree,
                                                                comm = get(paste0("community_1_bcr_", i, "_fil")))$phy)
  assign(paste0("passerine.picCleanComm.", i), match.phylo.comm(phy = passerine.tree,
                                                                comm = get(paste0("community_1_bcr_", i, "_fil")))$comm)
}

# create phylogenetic distance matrices
for(i in 1:37){
  if(nrow(get(paste0("wewpew_1_bcr_", i, "_df"))) > 0){
    assign(paste0("passerine.cophenDist.", i), cophenetic.phylo(get(paste0("passerine.picCleanTree.", i))))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}


#subsample dataframes
for(i in 1:37){
  input_df_name <- paste0("passerine.picCleanComm.", i)
  subset_name <- paste0("passerine.1.picCleanComm.", i, ".subsample")
  if(nrow(get(input_df_name)) > 0){
    #subset the dataframe using rownames
    subset_df <- get(input_df_name)[rownames(get(input_df_name)) %in% intsample, ]
    #assign the new subsetted dataframe
    assign(subset_name, subset_df)
  } else {
    message(paste0("Skipped BCR ", i)) # show skipped BCRs
  }
}

#run NRI on each BCR region subsample
for(i in 1:37){
  if(nrow(get(paste0("wewpew_1_bcr_", i, "_df"))) > 0){
    assign(paste0("passerine.ses.mpd.1.bcr", i), ses.mpd(get(paste0("passerine.1.picCleanComm.", i, ".subsample")),
                                                         get(paste0("passerine.cophenDist.", i)),
                                                         null.model = "taxa.labels",
                                                         abundance.weighted = FALSE,
                                                         runs = 100))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}

#
# #export tables for each BCR's NRI results
for(i in 1:37){
  if(nrow(get(paste0("wewpew_1_bcr_", i, "_df"))) > 0){
    write.csv(get(paste0("passerine.ses.mpd.1.bcr", i)), paste0("output/bcr_tables/1/passerine.ses.mpd.1.bcr.", i, ".redo.csv"))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}


### replicate process with NTI

for(i in 1:37){
  if(nrow(get(paste0("wewpew_1_bcr_", i, "_df"))) > 0){
    assign(paste0("passerine.ses.mntd.1.bcr", i), ses.mntd(get(paste0("passerine.1.picCleanComm.", i, ".subsample")),
                                                           get(paste0("passerine.cophenDist.", i)),
                                                           null.model = "taxa.labels",
                                                           abundance.weighted = FALSE,
                                                           runs = 100))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}


#export tables for each BCR's NTI results
for(i in 1:37){
  if(nrow(get(paste0("wewpew_1_bcr_", i, "_df"))) > 0){
    write.csv(get(paste0("passerine.ses.mntd.1.bcr", i)), paste0("output/bcr_tables/1/passerine.ses.mntd.1.bcr.", i, ".redo.csv"))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}


#
#re-import tables
for(i in 1:37){
  if(file.exists(paste0("output/bcr_tables/1/passerine.ses.mntd.1.bcr.", i, ".csv"))){
    assign(paste0("passerine.ses.mntd.1.bcr", i), read.csv(paste0("output/bcr_tables/1/passerine.ses.mntd.1.bcr.", i, ".redo.csv")))
    assign(paste0("passerine.ses.mpd.1.bcr", i), read.csv(paste0("output/bcr_tables/1/passerine.ses.mpd.1.bcr.", i, ".redo.csv")))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}

#create joint visualizations

#create BCR-specific row columns
for(i in 1:37){
  if(nrow(get(paste0("wewpew_1_bcr_", i, "_df"))) > 0){
    assign(paste0("passerine.ses.mpd.1.bcr", i), read.csv(paste0("output/bcr_tables/1/passerine.ses.mpd.1.bcr.", i, ".redo.csv"), row.names=1, header=TRUE))
    assign(paste0("passerine.ses.mpd.1.bcr", i, ".rows"), as.numeric(row.names(get(paste0("passerine.ses.mpd.1.bcr", i)))))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}

#create a list of all the NRI dataframes
week1bcrlist <- list(passerine.ses.mpd.1.bcr5, passerine.ses.mpd.1.bcr6,
                     passerine.ses.mpd.1.bcr7, passerine.ses.mpd.1.bcr8, passerine.ses.mpd.1.bcr9, passerine.ses.mpd.1.bcr10,
                     passerine.ses.mpd.1.bcr11, passerine.ses.mpd.1.bcr12, passerine.ses.mpd.1.bcr13, passerine.ses.mpd.1.bcr14,
                     passerine.ses.mpd.1.bcr15, passerine.ses.mpd.1.bcr16, passerine.ses.mpd.1.bcr17, passerine.ses.mpd.1.bcr18,
                     passerine.ses.mpd.1.bcr19, passerine.ses.mpd.1.bcr20, passerine.ses.mpd.1.bcr21, passerine.ses.mpd.1.bcr22,
                     passerine.ses.mpd.1.bcr23, passerine.ses.mpd.1.bcr24, passerine.ses.mpd.1.bcr25, passerine.ses.mpd.1.bcr26,
                     passerine.ses.mpd.1.bcr27, passerine.ses.mpd.1.bcr28, passerine.ses.mpd.1.bcr29, passerine.ses.mpd.1.bcr30,
                     passerine.ses.mpd.1.bcr31, passerine.ses.mpd.1.bcr32, passerine.ses.mpd.1.bcr33, passerine.ses.mpd.1.bcr34,
                     passerine.ses.mpd.1.bcr35, passerine.ses.mpd.1.bcr36, passerine.ses.mpd.1.bcr37)

#create list of all row name vectors
week1bcrlistrows <- list(passerine.ses.mpd.1.bcr5.rows, passerine.ses.mpd.1.bcr6.rows,
                         passerine.ses.mpd.1.bcr7.rows, passerine.ses.mpd.1.bcr8.rows, passerine.ses.mpd.1.bcr9.rows, passerine.ses.mpd.1.bcr10.rows,
                         passerine.ses.mpd.1.bcr11.rows, passerine.ses.mpd.1.bcr12.rows, passerine.ses.mpd.1.bcr13.rows, passerine.ses.mpd.1.bcr14.rows,
                         passerine.ses.mpd.1.bcr15.rows, passerine.ses.mpd.1.bcr16.rows, passerine.ses.mpd.1.bcr17.rows, passerine.ses.mpd.1.bcr18.rows,
                         passerine.ses.mpd.1.bcr19.rows, passerine.ses.mpd.1.bcr20.rows, passerine.ses.mpd.1.bcr21.rows, passerine.ses.mpd.1.bcr22.rows,
                         passerine.ses.mpd.1.bcr23.rows, passerine.ses.mpd.1.bcr24.rows, passerine.ses.mpd.1.bcr25.rows, passerine.ses.mpd.1.bcr26.rows,
                         passerine.ses.mpd.1.bcr27.rows, passerine.ses.mpd.1.bcr28.rows, passerine.ses.mpd.1.bcr29.rows, passerine.ses.mpd.1.bcr30.rows,
                         passerine.ses.mpd.1.bcr31.rows, passerine.ses.mpd.1.bcr32.rows, passerine.ses.mpd.1.bcr33.rows, passerine.ses.mpd.1.bcr34.rows,
                         passerine.ses.mpd.1.bcr35.rows, passerine.ses.mpd.1.bcr36.rows, passerine.ses.mpd.1.bcr37.rows)


#add columns for x coords, y coords, NRI data
for(i in 1:length(week1bcrlist)){
  week1bcrlist[[i]]$x <- rep(NA, nrow(week1bcrlist[[i]]))
  week1bcrlist[[i]]$y <- rep(NA, nrow(week1bcrlist[[i]]))
  week1bcrlist[[i]]$NRI <- rep(NA, nrow(week1bcrlist[[i]]))
  week1bcrlist[[i]]$BCR <- rep(NA, nrow(week1bcrlist[[i]]))
}

#
# #use agg spatraster of minimum length df
wewpew_1_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 1), us), crs, method = "near"), fact=2)

#add coords and NRI to each subsampled community
# #valid BCRs in week 1
validBCR1 <- c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,1,27,28,29,30,31,32,33,34,35,36,37)


for(i in 1:length(week1bcrlist)){
  for(j in 1:nrow(week1bcrlist[[i]])){
    week1bcrlist[[i]]$x[j] <- xyFromCell(wewpew_1_agg, week1bcrlistrows[[i]][j])[,1]
    week1bcrlist[[i]]$y[j] <- xyFromCell(wewpew_1_agg, week1bcrlistrows[[i]][j])[,2]
  }
  week1bcrlist[[i]]$NRI <- -1*week1bcrlist[[i]]$mpd.obs.z
  week1bcrlist[[i]]$BCR <- validBCR1[i]
  week1bcrlist[[i]] <- tibble::rownames_to_column(week1bcrlist[[i]], "community")
}


#combine into one dataframe
combined_week_1_nri <- bind_rows(week1bcrlist)

#plot community NRI
country_map <- ggplot()+
  geom_sf(data = countries) +
  coord_sf(xlim=c(-2750000, 3000000), ylim=c(-2000000, 4250000))


pdf(file = 'output/figures/1/passerine.1.BCR.NRI.pdf')
country_map +
  geom_point(data = combined_week_1_nri, alpha = 0.5, aes(x = x, y = y, color = NRI)) +
  scale_color_gradient2("NRI", low = "red", high = "blue4",
                        mid = "gray90",
                        midpoint = 0, limits = c(-4,5),
                        na.value = "gray97") +
  ggtitle("Jun 28") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "")
dev.off()
#
# #
# # repeat process with NTI
#
#create BCR-specific row columns
for(i in 1:37){
  if(nrow(get(paste0("wewpew_1_bcr_", i, "_df"))) > 0){
    assign(paste0("passerine.ses.mntd.1.bcr", i), read.csv(paste0("output/bcr_tables/1/passerine.ses.mntd.1.bcr.", i, ".redo.csv"), row.names=1, header=TRUE))
    assign(paste0("passerine.ses.mntd.1.bcr", i, ".rows"), as.numeric(row.names(get(paste0("passerine.ses.mntd.1.bcr", i)))))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}

#create a list of all the NTI dataframes
week1bcrlistNTI <- list(passerine.ses.mntd.1.bcr5, passerine.ses.mntd.1.bcr6,
                        passerine.ses.mntd.1.bcr7, passerine.ses.mntd.1.bcr8, passerine.ses.mntd.1.bcr9, passerine.ses.mntd.1.bcr10,
                        passerine.ses.mntd.1.bcr11, passerine.ses.mntd.1.bcr12, passerine.ses.mntd.1.bcr13, passerine.ses.mntd.1.bcr14,
                        passerine.ses.mntd.1.bcr15, passerine.ses.mntd.1.bcr16, passerine.ses.mntd.1.bcr17, passerine.ses.mntd.1.bcr18,
                        passerine.ses.mntd.1.bcr19, passerine.ses.mntd.1.bcr20, passerine.ses.mntd.1.bcr21, passerine.ses.mntd.1.bcr22,
                        passerine.ses.mntd.1.bcr23, passerine.ses.mntd.1.bcr24, passerine.ses.mntd.1.bcr25, passerine.ses.mntd.1.bcr26,
                        passerine.ses.mntd.1.bcr27, passerine.ses.mntd.1.bcr28, passerine.ses.mntd.1.bcr29, passerine.ses.mntd.1.bcr30,
                        passerine.ses.mntd.1.bcr31, passerine.ses.mntd.1.bcr32, passerine.ses.mntd.1.bcr33, passerine.ses.mntd.1.bcr34,
                        passerine.ses.mntd.1.bcr35, passerine.ses.mntd.1.bcr36, passerine.ses.mntd.1.bcr37)

#create list of all row name vectors
week1bcrlistrowsNTI <- list(passerine.ses.mntd.1.bcr5.rows, passerine.ses.mntd.1.bcr6.rows,
                            passerine.ses.mntd.1.bcr7.rows, passerine.ses.mntd.1.bcr8.rows, passerine.ses.mntd.1.bcr9.rows, passerine.ses.mntd.1.bcr10.rows,
                            passerine.ses.mntd.1.bcr11.rows, passerine.ses.mntd.1.bcr12.rows, passerine.ses.mntd.1.bcr13.rows, passerine.ses.mntd.1.bcr14.rows,
                            passerine.ses.mntd.1.bcr15.rows, passerine.ses.mntd.1.bcr16.rows, passerine.ses.mntd.1.bcr17.rows, passerine.ses.mntd.1.bcr18.rows,
                            passerine.ses.mntd.1.bcr19.rows, passerine.ses.mntd.1.bcr20.rows, passerine.ses.mntd.1.bcr21.rows, passerine.ses.mntd.1.bcr22.rows,
                            passerine.ses.mntd.1.bcr23.rows, passerine.ses.mntd.1.bcr24.rows, passerine.ses.mntd.1.bcr25.rows, passerine.ses.mntd.1.bcr26.rows,
                            passerine.ses.mntd.1.bcr27.rows, passerine.ses.mntd.1.bcr28.rows, passerine.ses.mntd.1.bcr29.rows, passerine.ses.mntd.1.bcr30.rows,
                            passerine.ses.mntd.1.bcr31.rows, passerine.ses.mntd.1.bcr32.rows, passerine.ses.mntd.1.bcr33.rows, passerine.ses.mntd.1.bcr34.rows,
                            passerine.ses.mntd.1.bcr35.rows, passerine.ses.mntd.1.bcr36.rows, passerine.ses.mntd.1.bcr37.rows)


#add columns for x coords, y coords, NTI data
for(i in 1:length(week1bcrlistNTI)){
  week1bcrlistNTI[[i]]$x <- rep(NA, nrow(week1bcrlistNTI[[i]]))
  week1bcrlistNTI[[i]]$y <- rep(NA, nrow(week1bcrlistNTI[[i]]))
  week1bcrlistNTI[[i]]$NTI <- rep(NA, nrow(week1bcrlistNTI[[i]]))
  week1bcrlistNTI[[i]]$BCR <- rep(NA, nrow(week1bcrlistNTI[[i]]))
}


# #use agg spatraster of minimum length df
wewpew_1_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 1), us), crs, method = "near"), fact=2)

#add coords and NTI to each subsampled community
for(i in 1:length(week1bcrlistNTI)){
  for(j in 1:nrow(week1bcrlistNTI[[i]])){
    week1bcrlistNTI[[i]]$x[j] <- xyFromCell(wewpew_1_agg, week1bcrlistrowsNTI[[i]][j])[,1]
    week1bcrlistNTI[[i]]$y[j] <- xyFromCell(wewpew_1_agg, week1bcrlistrowsNTI[[i]][j])[,2]
  }
  week1bcrlistNTI[[i]]$NTI <- -1*week1bcrlistNTI[[i]]$mntd.obs.z
  week1bcrlistNTI[[i]]$BCR <- validBCR1[i]
  week1bcrlistNTI[[i]] <- tibble::rownames_to_column(week1bcrlistNTI[[i]], "community")
}

#combine into one dataframe
combined_week_1_NTI <- bind_rows(week1bcrlistNTI)

#plot community NTI
country_map <- ggplot()+
  geom_sf(data = countries) +
  coord_sf(xlim=c(-2750000, 3000000), ylim=c(-2000000, 4250000))


pdf(file = 'output/figures/1/passerine.1.BCR.NTI.pdf')
country_map +
  geom_point(data = combined_week_1_NTI, alpha = 0.5, aes(x = x, y = y, color = NTI)) +
  scale_color_gradient2("NTI", low = "red", high = "blue4",
                        mid = "gray90",
                        midpoint = 0, limits = c(-5,5),
                        na.value = "gray97") +
  ggtitle("Jun 28") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "")
dev.off()
#
#

#re-create row names
for(i in 1:length(week1bcrlist)){
  assign(paste0("passerine.ses.mpd.1.bcr", validBCR1[i], ".rows"), week1bcrlist[[i]][,1])
  assign(paste0("passerine.ses.mntd.1.bcr", validBCR1[i], ".rows"), week1bcrlistNTI[[i]][,1])
}

#export combined NRI and NTI dataframes
write.csv(combined_week_1_nri, "output/bcr_tables/1/combined_NRI_1_redo.csv")
write.csv(combined_week_1_NTI, "output/bcr_tables/1/combined_NTI_1_redo.csv")

##### replicate process with continental pool tables
#re-import tables
passerine.ses.mntd.1.continental <- read.csv("output/tables/1/passerine.ses.mntd.1_redo.csv", row.names=1)
passerine.ses.mpd.1.continental <- read.csv("output/tables/1/passerine.ses.mpd.1_redo.csv", row.names=1)

passerine.ses.mpd.1.continental.rows <- as.numeric(row.names(passerine.ses.mpd.1.continental))

#NRI continental
#add columns for x coords, y coords, NRI data
passerine.ses.mpd.1.continental$x <- rep(NA, nrow(passerine.ses.mpd.1.continental))
passerine.ses.mpd.1.continental$y <- rep(NA, nrow(passerine.ses.mpd.1.continental))
passerine.ses.mpd.1.continental$NRI <- rep(NA, nrow(passerine.ses.mpd.1.continental))
passerine.ses.mpd.1.continental$BCR <- rep(NA, nrow(passerine.ses.mpd.1.continental))

# #use agg spatraster of minimum length df
wewpew_1_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 1), us), crs, method = "near"), fact=2)

#populate x and y
for(i in 1:nrow(passerine.ses.mpd.1.continental)){
  passerine.ses.mpd.1.continental$x[i] <- xyFromCell(wewpew_1_agg, passerine.ses.mpd.1.continental.rows[i])[,1]
  passerine.ses.mpd.1.continental$y[i] <- xyFromCell(wewpew_1_agg, passerine.ses.mpd.1.continental.rows[i])[,2]
}

passerine.ses.mpd.1.continental$NRI <- -1*passerine.ses.mpd.1.continental$mpd.obs.z
passerine.ses.mpd.1.continental <- tibble::rownames_to_column(passerine.ses.mpd.1.continental, "community")

#create list of BCR rows
bcrrowlist <- list(bcr1rows, bcr2rows, bcr3rows, bcr4rows, bcr5rows, bcr6rows, bcr7rows, bcr8rows,
                   bcr9rows, bcr10rows, bcr11rows, bcr12rows, bcr13rows, bcr14rows, bcr15rows, bcr16rows,
                   bcr17rows, bcr18rows, bcr19rows, bcr20rows, bcr21rows, bcr22rows, bcr23rows, bcr24rows,
                   bcr25rows, bcr26rows, bcr27rows, bcr28rows, bcr29rows, bcr30rows, bcr31rows, bcr32rows,
                   bcr33rows, bcr34rows, bcr35rows, bcr36rows, bcr37rows)

#match BCRs to continental community
for(i in 1:nrow(passerine.ses.mpd.1.continental)){
  if(any(sapply(bcrrowlist, function(x) passerine.ses.mpd.1.continental[i, "community"] %in% x)) == TRUE){
    passerine.ses.mpd.1.continental[i, "BCR"] <- which(sapply(bcrrowlist, function(x) passerine.ses.mpd.1.continental[i, "community"] %in% x))
  }
  else {
    print("Pass")
  }
}

migratory_proportion_df_week_1 <- read.csv("output/community_data/migratory_proportion_week_1.csv")

passerine.ses.mpd.1.continental$proportion <- rep(NA, nrow(passerine.ses.mpd.1.continental))
for(i in 1:nrow(passerine.ses.mpd.1.continental)){ 
  passerine.ses.mpd.1.continental[i,"proportion"] <- migratory_proportion_df_week_1[migratory_proportion_df_week_1$community == passerine.ses.mpd.1.continental[i,"community"],"proportion"]
}







# 
# #plot community NRI
# country_map <- ggplot()+
#   geom_sf(data = countries) +
#   coord_sf(xlim=c(-2750000, 3000000), ylim=c(-2000000, 4250000))
# 
# 
# pdf(file = 'output/figures/1/passerine.1.BCR.NRI.pdf')
# country_map + 
#   geom_point(data = combined_week_1_nri, alpha = 0.5, aes(x = x, y = y, color = NRI)) +
#   scale_color_gradient2("NRI", low = "red", high = "blue4",
#                         mid = "gray90",
#                         midpoint = 0, limits = c(-4,5),
#                         na.value = "gray97") +
#   ggtitle("Jan 4") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5)) +
#   labs(x = "", y = "")
# dev.off()
# 
# 
# 
# # repeat process with NTI
# 
# #create BCR-specific row columns
# for(i in 1:37){
#   if(nrow(get(paste0("wewpew_1_bcr_", i, "_df"))) > 0){
#     assign(paste0("passerine.ses.mntd.1.bcr", i), read.csv(paste0("output/bcr_tables/1/passerine.ses.mntd.1.bcr.", i, ".csv"), row.names=1, header=TRUE))
#     assign(paste0("passerine.ses.mntd.1.bcr", i, ".rows"), as.numeric(row.names(get(paste0("passerine.ses.mntd.1.bcr", i)))))
#   }
#   else {
#     message(paste0("Skipped BCR ", i)) #show skipped BCRs
#   }
# }
# 
# #create a list of all the NTI dataframes
# week1bcrlistNTI <- list(passerine.ses.mntd.1.bcr5, passerine.ses.mntd.1.bcr6,
#                         passerine.ses.mntd.1.bcr7, passerine.ses.mntd.1.bcr8, passerine.ses.mntd.1.bcr9, passerine.ses.mntd.1.bcr10,
#                         passerine.ses.mntd.1.bcr11, passerine.ses.mntd.1.bcr12, passerine.ses.mntd.1.bcr13, passerine.ses.mntd.1.bcr14,
#                         passerine.ses.mntd.1.bcr15, passerine.ses.mntd.1.bcr16, passerine.ses.mntd.1.bcr17, passerine.ses.mntd.1.bcr18,
#                         passerine.ses.mntd.1.bcr19, passerine.ses.mntd.1.bcr20, passerine.ses.mntd.1.bcr21, passerine.ses.mntd.1.bcr22,
#                         passerine.ses.mntd.1.bcr23, passerine.ses.mntd.1.bcr24, passerine.ses.mntd.1.bcr25, passerine.ses.mntd.1.bcr26,
#                         passerine.ses.mntd.1.bcr27, passerine.ses.mntd.1.bcr28, passerine.ses.mntd.1.bcr29, passerine.ses.mntd.1.bcr30,
#                         passerine.ses.mntd.1.bcr31, passerine.ses.mntd.1.bcr32, passerine.ses.mntd.1.bcr33, passerine.ses.mntd.1.bcr34,
#                         passerine.ses.mntd.1.bcr35, passerine.ses.mntd.1.bcr36, passerine.ses.mntd.1.bcr37)
# 
# #create list of all row name vectors
# week1bcrlistrowsNTI <- list(passerine.ses.mntd.1.bcr5.rows, passerine.ses.mntd.1.bcr6.rows,
#                             passerine.ses.mntd.1.bcr7.rows, passerine.ses.mntd.1.bcr8.rows, passerine.ses.mntd.1.bcr9.rows, passerine.ses.mntd.1.bcr10.rows,
#                             passerine.ses.mntd.1.bcr11.rows, passerine.ses.mntd.1.bcr12.rows, passerine.ses.mntd.1.bcr13.rows, passerine.ses.mntd.1.bcr14.rows,
#                             passerine.ses.mntd.1.bcr15.rows, passerine.ses.mntd.1.bcr16.rows, passerine.ses.mntd.1.bcr17.rows, passerine.ses.mntd.1.bcr18.rows,
#                             passerine.ses.mntd.1.bcr19.rows, passerine.ses.mntd.1.bcr20.rows, passerine.ses.mntd.1.bcr21.rows, passerine.ses.mntd.1.bcr22.rows,
#                             passerine.ses.mntd.1.bcr23.rows, passerine.ses.mntd.1.bcr24.rows, passerine.ses.mntd.1.bcr25.rows, passerine.ses.mntd.1.bcr26.rows,
#                             passerine.ses.mntd.1.bcr27.rows, passerine.ses.mntd.1.bcr28.rows, passerine.ses.mntd.1.bcr29.rows, passerine.ses.mntd.1.bcr30.rows,
#                             passerine.ses.mntd.1.bcr31.rows, passerine.ses.mntd.1.bcr32.rows, passerine.ses.mntd.1.bcr33.rows, passerine.ses.mntd.1.bcr34.rows,
#                             passerine.ses.mntd.1.bcr35.rows, passerine.ses.mntd.1.bcr36.rows, passerine.ses.mntd.1.bcr37.rows)
# 
# 
# #add columns for x coords, y coords, NTI data
# for(i in 1:length(week1bcrlistNTI)){
#   week1bcrlistNTI[[i]]$x <- rep(NA, nrow(week1bcrlistNTI[[i]]))
#   week1bcrlistNTI[[i]]$y <- rep(NA, nrow(week1bcrlistNTI[[i]]))
#   week1bcrlistNTI[[i]]$NTI <- rep(NA, nrow(week1bcrlistNTI[[i]]))
# }
# 
# 
# # #use agg spatraster of minimum length df
# wewpew_1_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 1), us), crs, method = "near"), fact=2)
# 
# #add coords and NTI to each subsampled community
# for(i in 1:length(week1bcrlistNTI)){
#   for(j in 1:nrow(week1bcrlistNTI[[i]])){
#     week1bcrlistNTI[[i]]$x[j] <- xyFromCell(wewpew_1_agg, week1bcrlistrowsNTI[[i]][j])[,1]
#     week1bcrlistNTI[[i]]$y[j] <- xyFromCell(wewpew_1_agg, week1bcrlistrowsNTI[[i]][j])[,2]
#   }
#   week1bcrlistNTI[[i]]$NTI <- -1*week1bcrlistNTI[[i]]$mntd.obs.z
# }
# 
# #combine into one dataframe
# combined_week_1_NTI <- bind_rows(week1bcrlistNTI, .id = "ntaxa")
# 
# #plot community NTI
# country_map <- ggplot()+
#   geom_sf(data = countries) +
#   coord_sf(xlim=c(-2750000, 3000000), ylim=c(-2000000, 4250000))
# 
# 
# pdf(file = 'output/figures/1/passerine.1.BCR.NTI.pdf')
# country_map + 
#   geom_point(data = combined_week_1_NTI, alpha = 0.5, aes(x = x, y = y, color = NTI)) +
#   scale_color_gradient2("NTI", low = "red", high = "blue4",
#                         mid = "gray90",
#                         midpoint = 0, limits = c(-5,5),
#                         na.value = "gray97") +
#   ggtitle("Jan 4") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5)) +
#   labs(x = "", y = "")
# dev.off()

# 
# #NTI visualization
# 
# passerine.ses.mntd.1.bcr12.rows <- as.numeric(row.names(passerine.ses.mntd.1.bcr12))
# 
# passerine.ses.mntd.1.bcr12$x <- rep(NA, nrow(passerine.ses.mntd.1.bcr12))
# passerine.ses.mntd.1.bcr12$y <- rep(NA, nrow(passerine.ses.mntd.1.bcr12))
# passerine.ses.mntd.1.bcr12$NTI <- rep(NA, nrow(passerine.ses.mntd.1.bcr12))
# 
# 
# #add coords to each subsampled community
# for(i in 1:nrow(passerine.ses.mntd.1.bcr12)){
#   passerine.ses.mntd.1.bcr12$x[i] <- xyFromCell(wewpew_1_bcr_12, passerine.ses.mntd.1.bcr12.rows[i])[,1]
#   passerine.ses.mntd.1.bcr12$y[i] <- xyFromCell(wewpew_1_bcr_12, passerine.ses.mntd.1.bcr12.rows[i])[,2]
# }
# 
# passerine.ses.mntd.1.bcr12$NTI <- -1*passerine.ses.mntd.1.bcr12$mntd.obs.z
# 
# 
# #plot community NTI
# ggplot(data = passerine.ses.mntd.1.bcr12, mapping = aes(x = x, y = y)) +
#   geom_point(alpha = 0.5, aes(color = NTI)) +
#   scale_color_gradient2("NTI", low = "red", high = "blue4",
#                         mid = "gray90",
#                         midpoint = 0, limits = c(min(passerine.ses.mntd.1.bcr12$NTI), max(passerine.ses.mntd.1.bcr12$NTI)),
#                         na.value = "gray97") +
#   ggtitle("Jan 4") +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         plot.title = element_text(hjust = 0.5)) +
#   labs(x = "", y = "")







# create 0608 community matrix
#community0608 <- data.frame(matrix(nrow=nrow(refwar_0608_agg_df),ncol=length(species_passerines)))
#rownames(community0608) <- rownames(refwar_0608_agg_df)
#colnames(community0608) <- species_passerines

#reduce to minimum shared length (represented by RFWA)
#for(i in 1:length(species_names)){
#  assign(paste0(species_names[i], "_0608_agg_df"), read.csv(paste0("output/species_data/0608/", species_names[i], "_0608_agg.csv"), row.names = 1))
#}

#for(i in 1:length(species_names)){
#  assign(paste0(species_names[i], "_0608_reduced"), 
#         get(paste0(species_names[i], "_0608_agg_df"))[rownames(get(paste0(species_names[i], "_0608_agg_df"))) %in% rownames(refwar_0608_agg_df),])
#  write.csv(get(paste0(species_names[i], "_0608_reduced")), file = paste0("output/species_data/0608/", species_names[i], "_0608_reduced.csv"))
#}

#read reduced dataframe
#for(i in 1:length(species_names)){
#  assign(paste0(species_names[i], "_0608_reduced"), read.csv(paste0("output/species_data/0608/", species_names[i], "_0608_reduced.csv"), row.names = 1))
#}




#community1 <- read.csv("output/community_data/community_1.csv")

# Using picante -- clean trees/comm
#passerine.picCleanTree <- match.phylo.comm(phy = passerine.tree, 
#                                           comm = community1)$phy
#passerine.picCleanComm <- match.phylo.comm(phy = passerine.tree, 
#                                           comm = community1)$comm

# Phylogenetic alpha diversity

#passerine.1.PD <- pd(samp = passerine.picCleanComm, 
#                        tree = passerine.picCleanTree,
#                        include.root = FALSE)
#write.csv(passerine.1.PD, "output/tables/passerine.1.PD.csv")
#passerine.1.PD <- read.csv("output/tables/passerine.1.PD.csv")


#### Randomizations and null models

#random subsample
#passerine.33.picCleanComm.subsample <- sample_n(passerine.picCleanComm, 25000) #feel free to modify
#passerine.33.ses.pd.subsample <- ses.pd(samp = passerine.33.picCleanComm.subsample, 
#                                       tree = passerine.picCleanTree, 
#                                       null.model = "taxa.labels", 
#                                       runs = 100)
#write.csv(passerine.33.ses.pd.subsample, "output/tables/33/passerine.pd.33.subsample.csv")
#passerine.1.ses.pd.subsample <- read.csv("output/tables/1/passerine.pd.1.subsample.csv", row.names=1, header=TRUE)
#head(passerine.1.ses.pd.subsample)

#hist(warbler.ses.pd.subsample$pd.obs.z)
#hist(warbler.ses.pd.subsample$ntaxa)
#plot(warbler.ses.pd.subsample$pd.obs.z ~ warbler.ses.pd.subsample$ntaxa)


# Create a Phylogenetic Distance Matrix
#passerine.cophenDist <- cophenetic.phylo(passerine.picCleanTree)


#NRI on summer community
#passerine.ses.mpd.33 <- ses.mpd(passerine.33.picCleanComm.subsample, 
#                               passerine.cophenDist, 
#                               null.model = "taxa.labels", 
#                               abundance.weighted = FALSE, 
#                               runs = 100)
#hist(passerine.ses.mpd.1$mpd.obs.z)
#hist(passerine.ses.mpd.1$ntaxa)
#plot(passerine.ses.mpd.1$mpd.obs.z ~ passerine.ses.mpd.1$ntaxa)

#write.csv(passerine.ses.mpd.33, "output/tables/33/passerine.ses.mpd.33.csv")
#passerine.ses.mpd.1 <- read.csv("output/tables/1/passerine.ses.mpd.1.csv", row.names=1, header=TRUE)



### Nearest Taxon Index Calculations

#NTI on summer
#passerine.ses.mntd.33 <- ses.mntd(passerine.33.picCleanComm.subsample, 
#                                 passerine.cophenDist, 
#                                 null.model = "taxa.labels", 
#                                 abundance.weighted = FALSE, 
#                                 runs = 100) 
#hist(passerine.ses.mntd.1$mntd.obs.z)
#plot(passerine.ses.mntd.1$mntd.obs.z ~ passerine.ses.mntd.1$ntaxa)
#write.csv(passerine.ses.mntd.33, "output/tables/33/passerine.ses.mntd.33.csv")

# visualizing for initial results

#hist(passerine.0608.ses.pd.subsample$pd.obs.z)
#taxarange <- c(1:89)
#meanbytaxa <- rep(NA, length(taxarange))
#mean(passerine.0608.ses.pd.subsample$pd.obs.z[passerine.0608.ses.pd.subsample$ntaxa == 20])
#for(i in 1:length(taxarange)){
#  meanbytaxa[i] <- mean(passerine.0608.ses.pd.subsample$pd.obs.z[passerine.0608.ses.pd.subsample$ntaxa == i])
#}
#plot(meanbytaxa ~ c(1:89))
#hist(passerine.0608.ses.pd.subsample$ntaxa)

#visualizing PD ~ SR relationship
#ggplot(data = passerine.0608.ses.pd.subsample[passerine.0608.ses.pd.subsample$ntaxa > 1,], aes(x = ntaxa, y = pd.obs.z)) + 
#  geom_point(aes(col = pd.obs.z), alpha = 0.05) +
#  scale_color_gradient2("PD", low = "red", high = "blue4", 
#                        mid = "gray90",
#                        midpoint = 0, limits = c(min(passerine.0608.ses.pd.subsample$pd.obs.z), max(passerine.0608.ses.pd.subsample$pd.obs.z)), 
#                        na.value = "gray97") +
#  labs(x = "Species Richness", y = "SES of Phylogenetic Diversity") +
#  xlim(0, max(passerine.0608.ses.pd.subsample$ntaxa)) +
#  theme(axis.title.y = element_text(margin = margin(r = 15)))

#passerine.0608.ses.pd.subsample.rows <- as.numeric(row.names(passerine.0608.ses.pd.subsample))
#head(passerine.0608.ses.pd.subsample)
#head(community0608)

#passerine.0608.ses.pd.subsample$x <- rep(NA, nrow(passerine.0608.ses.pd.subsample))
#passerine.0608.ses.pd.subsample$y <- rep(NA, nrow(passerine.0608.ses.pd.subsample))

#for(i in 1:nrow(passerine.0608.ses.pd.subsample)){
#  passerine.0608.ses.pd.subsample$x[i] <- xyFromCell(refwar_0608_agg, passerine.0608.ses.pd.subsample.rows[i])[,1]
#  passerine.0608.ses.pd.subsample$y[i] <- xyFromCell(refwar_0608_agg, passerine.0608.ses.pd.subsample.rows[i])[,2]
#}




#plot(passerine.0608.ses.pd.subsample$y ~ passerine.0608.ses.pd.subsample$x)

# visualize PD in 0608 community

#map of PD
#pdf(file = "output/figures/warbler.0608.PD.pdf")
#ggplot(data = passerine.0608.ses.pd.subsample, mapping = aes(x = x, y = y)) + 
#  geom_point(alpha = 0.5, aes(color = pd.obs.z)) + 
#  scale_color_gradient2("SES PD", low = "red", high = "blue4", 
#                        mid = "gray90",
#                        midpoint = 0, limits = c(min(passerine.0608.ses.pd.subsample$pd.obs.z), max(passerine.0608.ses.pd.subsample$pd.obs.z)), 
#                        na.value = "gray97") +  
#  theme(axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank(),
#        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        plot.title = element_text(hjust = 0.5)) +
#  labs(x = "", y = "")
#dev.off()
