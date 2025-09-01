#Manuscript figures

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
#devtools::install_github("ropensci/rnaturalearthhires")
extract <- terra::extract
project <- terra::project #possible source of conf <- ion regarding rgdal for later

# Pure intro stuff

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

# download list of passerines

species_codes <- c("rotbec", "nobtyr", "olsfly", "grepew", "wewpew", "eawpew", "yebfly", "acafly",
                   "aldfly", "wilfly", "leafly", "hamfly", "gryfly", "dusfly", "pasfly", "corfly", 
                   "bubfly", "blkpho", "easpho", "saypho", "verfly", "ducfly", "astfly", "grcfly", 
                   "bncfly", "grekis", "subfly", "trokin", "coukin", "caskin", "thbkin", "weskin", 
                   "easkin", "grykin", "sctfly", "bkcvir1", "whevir", "belvir", "gryvir", "hutvir", 
                   "yetvir", "casvir", "buhvir", "plsvir", "phivir", "warvir", "reevir1", "bkwvir", 
                   "logshr", "norshr4", "gryjay", "brnjay", "grnjay", "pinjay", "stejay", "blujay", 
                   "flsjay", "issjay", "cowscj1", "wooscj2", "mexjay4", "bkbmag1", "yebmag", "clanut", 
                   "amecro", "fiscro", "chirav", "comrav", "carchi", "bkcchi", "mouchi", "mexchi", 
                   "chbchi", "borchi2", "gyhchi", "britit", "oaktit", "juntit1", "tuftit", "blctit4", 
                   "verdin", "horlar", "nrwswa", "purmar", "treswa", "vigswa", "banswa", "cliswa", 
                   "cavswa", "rewbul", "arcwar1", "bushti", "wrenti", "ruckin", "gockin", "rebnut", 
                   "whbnut", "pygnut", "bnhnut", "brncre", "buggna", "bktgna", "calgna", "bkcgna", 
                   "rocwre", "canwre", "houwre", "pacwre1", "winwre3", "sedwre1", "marwre", "carwre", 
                   "bewwre", "cacwre", "amedip", "eursta", "commyn", "grycat", "cubthr", "brnthr", 
                   "lobthr", "benthr", "calthr", "lecthr", "crithr", "sagthr", "normoc", "easblu", 
                   "wesblu", "moublu", "towsol", "varthr", "veery", "gycthr", "bicthr", "swathr", 
                   "herthr", "woothr", "clcrob", "amerob", "norwhe", "bohwax", "cedwax", "phaino", 
                   "oliwar", "orabis1", "nutman", "houspa", "eutspa", "eaywag", "amepip", "sprpip", 
                   "evegro", "pingro", "gcrfin", "bkrfin", "bcrfin", "houfin", "purfin", "casfin", 
                   "comred", "hoared", "redcro", "whwcro", "pinsis", "lesgol", "lawgol", "amegfi", 
                   "laplon", "chclon", "smilon", "mcclon", "snobun", "mckbun", "ruwspa", "botspa", 
                   "casspa", "bacspa", "graspa", "olispa", "chispa", "clcspa", "bkcspa", "fiespa", 
                   "brespa", "fisspa", "bktspa", "larspa", "larbun", "amtspa", "foxspa", "daejun", 
                   "yeejun", "whcspa", "gocspa", "harspa", "whtspa", "sagspa1", "belspa2", "vesspa", 
                   "lecspa", "seaspa", "nstspa", "sstspa", "savspa", "baispa", "henspa", "sonspa", 
                   "linspa", "swaspa", "cantow", "abetow", "caltow", "rucspa", "gnttow", "spotow", 
                   "eastow", "yebcha", "yehbla", "boboli", "wesmea", "easmea", "orcori", "hooori", 
                   "bulori", "spbori", "altori", "audori", "balori", "scoori", "rewbla", "tribla", 
                   "shicow", "brocow", "bnhcow", "rusbla", "brebla", "comgra", "botgra", "grtgra", 
                   "ovenbi1", "woewar1", "louwat", "norwat", "gowwar", "buwwar", "bawwar", "prowar",
                   "swawar", "tenwar", "orcwar", "colwar", "lucwar", "naswar", "virwar", "conwar", 
                   "macwar", "mouwar", "kenwar", "comyel", "hoowar", "amered", "kirwar", "camwar", 
                   "cerwar", "norpar", "tropar", "magwar", "babwar", "bkbwar", "yelwar", "chswar", 
                   "bkpwar", "btbwar", "palwar", "pinwar", "yerwar", "yetwar", "prawar", "grawar", 
                   "btywar", "towwar", "herwar", "gchwar", "btnwar", "canwar", "wlswar", "refwar", 
                   "paired", "heptan", "sumtan", "scatan", "westan", "norcar", "pyrrhu", "robgro", 
                   "bkhgro", "blugrb1", "lazbun", "indbun", "varbun", "paibun", "dickci")

# "whcsee1", "blueth

species_names <- c("rotbec", "nobtyr", "olsfly", "grepew", "wewpew", "eawpew", "yebfly", "acafly",
                   "aldfly", "wilfly", "leafly", "hamfly", "gryfly", "dusfly", "pasfly", "corfly", 
                   "bubfly", "blkpho", "easpho", "saypho", "verfly", "ducfly", "astfly", "grcfly", 
                   "bncfly", "grekis", "subfly", "trokin", "coukin", "caskin", "thbkin", "weskin", 
                   "easkin", "grykin", "sctfly", "bkcvir1", "whevir", "belvir", "gryvir", "hutvir", 
                   "yetvir", "casvir", "buhvir", "plsvir", "phivir", "warvir", "reevir1", "bkwvir", 
                   "logshr", "norshr4", "gryjay", "brnjay", "grnjay", "pinjay", "stejay", "blujay", 
                   "flsjay", "issjay", "cowscj1", "wooscj2", "mexjay4", "bkbmag1", "yebmag", "clanut", 
                   "amecro", "fiscro", "chirav", "comrav", "carchi", "bkcchi", "mouchi", "mexchi", 
                   "chbchi", "borchi2", "gyhchi", "britit", "oaktit", "juntit1", "tuftit", "blctit4", 
                   "verdin", "horlar", "nrwswa", "purmar", "treswa", "vigswa", "banswa", "cliswa", 
                   "cavswa", "rewbul", "arcwar1", "bushti", "wrenti", "ruckin", "gockin", "rebnut", 
                   "whbnut", "pygnut", "bnhnut", "brncre", "buggna", "bktgna", "calgna", "bkcgna", 
                   "rocwre", "canwre", "houwre", "pacwre1", "winwre3", "sedwre1", "marwre", "carwre", 
                   "bewwre", "cacwre", "amedip", "eursta", "commyn", "grycat", "cubthr", "brnthr", 
                   "lobthr", "benthr", "calthr", "lecthr", "crithr", "sagthr", "normoc", "easblu", 
                   "wesblu", "moublu", "towsol", "varthr", "veery", "gycthr", "bicthr", "swathr", 
                   "herthr", "woothr", "clcrob", "amerob", "norwhe", "bohwax", "cedwax", "phaino", 
                   "oliwar", "orabis1", "nutman", "houspa", "eutspa", "eaywag", "amepip", "sprpip", 
                   "evegro", "pingro", "gcrfin", "bkrfin", "bcrfin", "houfin", "purfin", "casfin", 
                   "comred", "hoared", "redcro", "whwcro", "pinsis", "lesgol", "lawgol", "amegfi", 
                   "laplon", "chclon", "smilon", "mcclon", "snobun", "mckbun", "ruwspa", "botspa", 
                   "casspa", "bacspa", "graspa", "olispa", "chispa", "clcspa", "bkcspa", "fiespa", 
                   "brespa", "fisspa", "bktspa", "larspa", "larbun", "amtspa", "foxspa", "daejun", 
                   "yeejun", "whcspa", "gocspa", "harspa", "whtspa", "sagspa1", "belspa2", "vesspa", 
                   "lecspa", "seaspa", "nstspa", "sstspa", "savspa", "baispa", "henspa", "sonspa", 
                   "linspa", "swaspa", "cantow", "abetow", "caltow", "rucspa", "gnttow", "spotow", 
                   "eastow", "yebcha", "yehbla", "boboli", "wesmea", "easmea", "orcori", "hooori", 
                   "bulori", "spbori", "altori", "audori", "balori", "scoori", "rewbla", "tribla", 
                   "shicow", "brocow", "bnhcow", "rusbla", "brebla", "comgra", "botgra", "grtgra", 
                   "ovenbi1", "woewar1", "louwat", "norwat", "gowwar", "buwwar", "bawwar", "prowar",
                   "swawar", "tenwar", "orcwar", "colwar", "lucwar", "naswar", "virwar", "conwar", 
                   "macwar", "mouwar", "kenwar", "comyel", "hoowar", "amered", "kirwar", "camwar", 
                   "cerwar", "norpar", "tropar", "magwar", "babwar", "bkbwar", "yelwar", "chswar", 
                   "bkpwar", "btbwar", "palwar", "pinwar", "yerwar", "yetwar", "prawar", "grawar", 
                   "btywar", "towwar", "herwar", "gchwar", "btnwar", "canwar", "wlswar", "refwar", 
                   "paired", "heptan", "sumtan", "scatan", "westan", "norcar", "pyrrhu", "robgro", 
                   "bkhgro", "blugrb1", "lazbun", "indbun", "varbun", "paibun", "dickci")

#, "whcsee1"

#download data for each species

#for(i in 1:length(species_codes)){
#  str_sub(species_codes[i]) <- ebirdst_download(species = species_codes[i])
#}

#create passerine abd seasonal mean files

#for(i in 1:length(species_names)){
#  assign(paste0(species_names[i], "_path"), get_species_path(species_names[i]))
#}

#load_raster(paste0("data/ebirdst/2021/", species_names[1]), product = "abundance",
#            period = "seasonal", metric = "mean", 
#            resolution = "hr")

#for(i in 1:length(species_names)){
#  assign(paste0(species_names[i], "_abd_seasonal_mean"), load_raster(paste0("data/ebirdst/2021/", species_names[i]), product = "abundance", 
#                                                                     period = "seasonal", metric = "mean", 
#                                                                     resolution = "hr"))
#}

#see if you can get a tree working

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

# Load all required libraries at once
lapply(required.libraries, require, character.only = TRUE)

#load in the tree
tree.passerine <- read.nexus("phylogeny_data/jetz_passerine_tree/jetz_consensus_tree.tre")
tree.passerine <- root(tree.passerine, outgroup = c(263:297))


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
#, "Sporophila_torqueola", "Luscinia svecica"

passerine.tree <- drop.tip(tree.passerine,tree.passerine$tip.label[-match(species_passerines, tree.passerine$tip.label)])

#now try to get to individual weeks for each warbler (wrap)
#this gets me to a specific week (0608)


ovenbi1_abd_seasonal_mean <- load_raster("data/ebirdst/2021/ovenbi1", product = "abundance", 
                                         period = "seasonal", metric = "mean", 
                                         resolution = "hr")

us <- ne_countries(country = "united states of america", type = "countries", returnclass = "sf") %>% 
  st_transform(st_crs(ovenbi1_abd_seasonal_mean))


### visualizations
#week 1
migratory_proportion_df_week_1 <- read.csv("output/community_data/migratory_proportion_week_1.csv")
migratory_proportion_df_week_1 <- data.frame(migratory_proportion_df_week_1[, c("community", "nmigrant", "nspecies", "proportion")],
                                             row.names = migratory_proportion_df_week_1$community) 
#use agg spatraster of minimum length df
wewpew_1_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/wewpew", product = "abundance", resolution = "hr"), 1), us), crs, method = "near"), fact=2)

migratory_proportion_df_week_1.rows <- as.numeric(migratory_proportion_df_week_1$community)

migratory_proportion_df_week_1$x <- rep(NA, nrow(migratory_proportion_df_week_1))
migratory_proportion_df_week_1$y <- rep(NA, nrow(migratory_proportion_df_week_1))

for(i in 1:nrow(migratory_proportion_df_week_1)){
  migratory_proportion_df_week_1$x[i] <- xyFromCell(wewpew_1_agg, migratory_proportion_df_week_1.rows[i])[,1]
  migratory_proportion_df_week_1$y[i] <- xyFromCell(wewpew_1_agg, migratory_proportion_df_week_1.rows[i])[,2]
}



country_map <- ggplot()+
  geom_sf(data = countries) +
  coord_sf(xlim=c(-2750000, 3000000), ylim=c(-2000000, 4250000))

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


#week 26
migratory_proportion_df_week_26 <- read.csv("output/community_data/migratory_proportion_week_26.csv")
migratory_proportion_df_week_26 <- data.frame(migratory_proportion_df_week_26[, c("community", "nmigrant", "nspecies", "proportion")],
                                              row.names = migratory_proportion_df_week_26$community) 
#use agg spatraster of minimum length df
eawpew_26_agg <- aggregate(project(crop(subset(load_raster("data/ebirdst/2021/eawpew", product = "abundance", resolution = "hr"), 26), us), crs, method = "near"), fact=2)

migratory_proportion_df_week_26.rows <- as.numeric(migratory_proportion_df_week_26$community)

migratory_proportion_df_week_26$x <- rep(NA, nrow(migratory_proportion_df_week_26))
migratory_proportion_df_week_26$y <- rep(NA, nrow(migratory_proportion_df_week_26))

for(i in 1:nrow(migratory_proportion_df_week_26)){
  migratory_proportion_df_week_26$x[i] <- xyFromCell(eawpew_26_agg, migratory_proportion_df_week_26.rows[i])[,1]
  migratory_proportion_df_week_26$y[i] <- xyFromCell(eawpew_26_agg, migratory_proportion_df_week_26.rows[i])[,2]
}



country_map <- ggplot()+
  geom_sf(data = countries) +
  coord_sf(xlim=c(-2750000, 3000000), ylim=c(-2000000, 4250000))

#pdf(file = 'output/figures/community_migratory_proportion_week_26.pdf')
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
#dev.off()


### sub in species richness

#week 1 species richness
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

#week 26 species richness
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

###histograms

#mig. prop. histogram
proportion_hist <- data.frame(matrix(NA, nrow = nrow(migratory_proportion_df_week_26)+nrow(migratory_proportion_df_week_1), ncol = 2))
colnames(proportion_hist) <- c("value", "type")
proportion_hist[c(1:nrow(migratory_proportion_df_week_1)), 1] <- migratory_proportion_df_week_1$proportion
proportion_hist[c((nrow(migratory_proportion_df_week_1) + 1)):(nrow(proportion_hist)),1] <- migratory_proportion_df_week_26$proportion
proportion_hist[c(1:nrow(migratory_proportion_df_week_1)), 2] <- "Non-Breeding"
proportion_hist[c((nrow(migratory_proportion_df_week_1) + 1)):(nrow(proportion_hist)),2] <- "Breeding"
tail(proportion_hist)

ggplot(proportion_hist, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("Proportion Species Migratory") +
  ylab("Density") + 
  guides(fill=guide_legend(title="Season")) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25))

#richness histogram
nspecies_hist <- data.frame(matrix(NA, nrow = nrow(migratory_proportion_df_week_26)+nrow(migratory_proportion_df_week_1), ncol = 2))
colnames(nspecies_hist) <- c("value", "type")
nspecies_hist[c(1:nrow(migratory_proportion_df_week_1)), 1] <- migratory_proportion_df_week_1$nspecies
nspecies_hist[c((nrow(migratory_proportion_df_week_1) + 1)):(nrow(nspecies_hist)),1] <- migratory_proportion_df_week_26$nspecies
nspecies_hist[c(1:nrow(migratory_proportion_df_week_1)), 2] <- "Non-Breeding"
nspecies_hist[c((nrow(migratory_proportion_df_week_1) + 1)):(nrow(nspecies_hist)),2] <- "Breeding"
tail(nspecies_hist)

ggplot(nspecies_hist, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("Species Richness") +
  ylab("Density") + 
  guides(fill=guide_legend(title="Season")) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25))


### get extreme values
#(load all_BCR values from annual cycle analyses script)

all_combined_week_NRIs_df <- data.table::rbindlist(all_combined_week_NRIs, use.names=TRUE, idcol=TRUE, fill=TRUE) %>% 
  as.data.frame()


all_combined_week_NTIs_df <- data.table::rbindlist(all_combined_week_NTIs, use.names=TRUE, idcol=TRUE, fill=TRUE) %>% 
  as.data.frame()

max(na.omit(all_combined_week_NRIs_df$NRI)) #7.109175
min(na.omit(all_combined_week_NRIs_df$NRI)) #-4.155872
median(na.omit(all_combined_week_NRIs_df$NRI)) #-0.2703786
sd(na.omit(all_combined_week_NRIs_df$NRI)) #1.039971
length(na.omit(all_combined_week_NRIs_df$NRI[all_combined_week_NRIs_df$NRI > 0]))/length(na.omit(all_combined_week_NRIs_df$NRI)) #0.389

max(na.omit(all_combined_week_NTIs_df$NTI)) #4.776999
min(na.omit(all_combined_week_NTIs_df$NTI)) #-4.575684
median(na.omit(all_combined_week_NTIs_df$NTI)) #0.3231434
sd(na.omit(all_combined_week_NTIs_df$NTI)) #0.973363
length(na.omit(all_combined_week_NTIs_df$NTI[all_combined_week_NTIs_df$NTI > 0]))/length(na.omit(all_combined_week_NTIs_df$NTI)) #0.626


max(na.omit(all_combined_week_NTIs_df$ntaxa)) #109
min(na.omit(all_combined_week_NTIs_df$ntaxa)) #0
median(na.omit(all_combined_week_NTIs_df$ntaxa)) #34
sd(na.omit(all_combined_week_NTIs_df$ntaxa)) #16.13621

max(na.omit(all_combined_week_NTIs_df$ntaxa))

#breeding
max(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26,"ntaxa"])) #84
min(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26,"ntaxa"])) #1
mean(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26,"ntaxa"])) #40.08
sd(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 26,"ntaxa"])) #11.62
mean(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 26,"proportion"])) #0.885
sd(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 26,"proportion"])) #0.076

#wintering
max(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1,"ntaxa"])) #78
min(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1,"ntaxa"])) #1
mean(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1,"ntaxa"])) #23.34
sd(na.omit(all_combined_week_NTIs_df[all_combined_week_NTIs_df$.id == 1,"ntaxa"])) #13.85
mean(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 1,"proportion"])) #0.748
sd(na.omit(all_combined_week_NRIs_df[all_combined_week_NRIs_df$.id == 1,"proportion"])) #0.151

#different BCRs
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

#### try to get nmigrant data again

for(i in 1:52){
  assign(paste0("week_", i, "_mig_prop"), read.csv(paste0("output/community_data/migratory_proportion_week_", i, ".csv"), row.names=1))
}

mig_prop_list <- list(week_1_mig_prop, week_2_mig_prop, week_3_mig_prop, week_4_mig_prop, week_5_mig_prop, week_6_mig_prop, week_7_mig_prop, week_8_mig_prop,
                week_9_mig_prop, week_10_mig_prop, week_11_mig_prop, week_12_mig_prop, week_13_mig_prop, week_14_mig_prop, week_15_mig_prop, week_16_mig_prop,
                week_17_mig_prop, week_18_mig_prop, week_19_mig_prop, week_20_mig_prop, week_21_mig_prop, week_22_mig_prop, week_23_mig_prop, week_24_mig_prop,
                week_25_mig_prop, week_26_mig_prop, week_27_mig_prop, week_28_mig_prop, week_29_mig_prop, week_30_mig_prop, week_31_mig_prop, week_32_mig_prop,
                week_33_mig_prop, week_34_mig_prop, week_35_mig_prop, week_36_mig_prop, week_37_mig_prop, week_38_mig_prop, week_39_mig_prop, week_40_mig_prop,
                week_41_mig_prop, week_42_mig_prop, week_43_mig_prop, week_44_mig_prop, week_45_mig_prop, week_46_mig_prop, week_47_mig_prop, week_48_mig_prop,
                week_49_mig_prop, week_50_mig_prop, week_51_mig_prop, week_52_mig_prop)

all_mig_prop_df <- data.table::rbindlist(mig_prop_list, use.names=TRUE, idcol=TRUE, fill=TRUE) %>% 
  as.data.frame()

#this dataframe can do statistics for first results paragraph

median(na.omit(all_mig_prop_df$proportion)) #0.857



#histogram for NRI and NTI values (BCR-level)
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




lmer_test <- lmer(NRI ~ ntaxa + latitude + (1 | BCR), data = all_combined_week_NRIs_df)
ranova(lmer_test)
tab_model(lmer_test)

lmer_test_NTI <- lmer(NTI ~ ntaxa + latitude + (1 | BCR), data = all_combined_week_NTIs_df)
ranova(lmer_test_NTI)
tab_model(lmer_test_NTI)


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
