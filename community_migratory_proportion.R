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


#### try to assign proportion migratory for each community

#load in the community
#load in minimum length df
wewpew_1_agg_df <- read.csv("output/species_data/1/wewpew_1_agg.csv", row.names = 1)
#read existing community table
#community1 <- data.table::fread("output/community_data/community_1.csv")
#assign rownames
#rownames(community1) <- rownames(wewpew_1_agg_df)
#rewrite community table with row names
#data.table::fwrite(community1, file = "output/community_data/community_1.csv", row.names = TRUE)
#reread community table
community1 <- data.table::fread("output/community_data/community_1.csv")
#assign correct rownames
community1 <- data.frame(community1[,-c("V1")], row.names = community1$V1) #you have to double check the number of rows here
head(community1)



#read in migratory status
migratory_status <- read.csv("migratory_status.csv")
#head(migratory_status)
#migratory_status[migratory_status$jetz_tree_latin_name == colnames(community1)[305],"migratory_status"] #returns migratory status for each species

#create empty migratory proportion dataframe
migratory_proportion_df_week_1 <- data.frame(community = rep(NA, nrow(community1)),
                                             nmigrant = rep(NA, nrow(community1)),
                                             nspecies = rep(NA, nrow(community1)),
                                             proportion = rep(NA, nrow(community1)))

#fill community numbers
migratory_proportion_df_week_1[,"community"] <- rownames(community1)

#fill diversity
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

#test example
#number_of_migrants(community1[224,] %>% select_if(~any(. > 0)))

#fill into dataset (will take a long time)
for(i in 1:nrow(migratory_proportion_df_week_1)){
  migratory_proportion_df_week_1[i,"nmigrant"] <- number_of_migrants(community1[i,] %>% select_if(~any(. > 0)))
}

migratory_proportion_df_week_1$proportion <- migratory_proportion_df_week_1$nmigrant/migratory_proportion_df_week_1$nspecies

write.csv(migratory_proportion_df_week_1, "output/community_data/migratory_proportion_week_1.csv")

#testrowcomm1 <- community1[1,] %>% select_if(~any(. > 0))

#colnames(community1[1,] %>% select_if(~any(. > 0))) #returns each indexed present species

#migratory_status[migratory_status$jetz_tree_latin_name == colnames(community1[1,] %>% select_if(~any(. > 0)))[1],"migratory_status"] #returns migratory status of each species


#list migratory species
#community1[1,][rowSums(community1[1,] > 0) == ncol(df), ]

#load in migratory proportion dataframes 
migratory_proportion_df_week_1 <- read.csv("output/community_data/migratory_proportion_week_1.csv")
migratory_proportion_df_week_1 <- data.frame(migratory_proportion_df_week_1[, c("community", "nmigrant", "nspecies", "proportion")],
                                             row.names = migratory_proportion_df_week_1$community) 


#extract by BCR

#week 1
for(i in 1:37){
  assign(paste0("bcr", i, "_mig_prop_week_1"), subset(migratory_proportion_df_week_1, row.names(migratory_proportion_df_week_1) %in% 
                                                        get(paste0("bcr", i, "rows"))))
}

#export by BCR
for(i in 1:37){
  if(nrow(get(paste0("bcr", i, "_mig_prop_week_1"))) > 0){
    write.csv(get(paste0("bcr", i, "_mig_prop_week_1")), paste0("output/bcr_community_data/1/bcr", i, "_mig_prop_week_1.csv"))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}

#week 2
migratory_proportion_df_week_2 <- read.csv("output/community_data/migratory_proportion_week_2.csv")
migratory_proportion_df_week_2 <- data.frame(migratory_proportion_df_week_2[, c("community", "nmigrant", "nspecies", "proportion")],
                                             row.names = migratory_proportion_df_week_2$community) 


for(i in 1:37){
  assign(paste0("bcr", i, "_mig_prop_week_2"), subset(migratory_proportion_df_week_2, row.names(migratory_proportion_df_week_2) %in% 
                                                        get(paste0("bcr", i, "rows"))))
}

#export by BCR
for(i in 1:37){
  if(nrow(get(paste0("bcr", i, "_mig_prop_week_2"))) > 0){
    write.csv(get(paste0("bcr", i, "_mig_prop_week_2")), paste0("output/bcr_community_data/2/bcr", i, "_mig_prop_week_2.csv"))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}


#week 16
migratory_proportion_df_week_16 <- read.csv("output/community_data/migratory_proportion_week_16.csv")
migratory_proportion_df_week_16 <- data.frame(migratory_proportion_df_week_16[, c("community", "nmigrant", "nspecies", "proportion")],
                                              row.names = migratory_proportion_df_week_16$community) 


for(i in 1:37){
  assign(paste0("bcr", i, "_mig_prop_week_16"), subset(migratory_proportion_df_week_16, row.names(migratory_proportion_df_week_16) %in% 
                                                         get(paste0("bcr", i, "rows"))))
}

#export by BCR
for(i in 1:37){
  if(nrow(get(paste0("bcr", i, "_mig_prop_week_16"))) > 0){
    write.csv(get(paste0("bcr", i, "_mig_prop_week_16")), paste0("output/bcr_community_data/16/bcr", i, "_mig_prop_week_16.csv"))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}


#week 19
migratory_proportion_df_week_19 <- read.csv("output/community_data/migratory_proportion_week_19.csv")
migratory_proportion_df_week_19 <- data.frame(migratory_proportion_df_week_19[, c("community", "nmigrant", "nspecies", "proportion")],
                                              row.names = migratory_proportion_df_week_19$community) 


for(i in 1:37){
  assign(paste0("bcr", i, "_mig_prop_week_19"), subset(migratory_proportion_df_week_19, row.names(migratory_proportion_df_week_19) %in% 
                                                         get(paste0("bcr", i, "rows"))))
}

#export by BCR
for(i in 1:37){
  if(nrow(get(paste0("bcr", i, "_mig_prop_week_19"))) > 0){
    write.csv(get(paste0("bcr", i, "_mig_prop_week_19")), paste0("output/bcr_community_data/19/bcr", i, "_mig_prop_week_19.csv"))
  }
  else {
    message(paste0("Skipped BCR ", i)) #show skipped BCRs
  }
}




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
        plot.title = element_text(hjust = 0.5)) +
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
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "")
dev.off()


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
        plot.title = element_text(hjust = 0.5)) +
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
        plot.title = element_text(hjust = 0.5)) +
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
  geom_density(alpha=0.5) +
  xlab("Proportion Species Migratory") +
  ylab("Density") + 
  guides(fill=guide_legend(title="Season")) +
  theme(legend.position = "bottom", text = element_text(size = 20))

#richness histogram
nspecies_hist <- data.frame(matrix(NA, nrow = nrow(migratory_proportion_df_week_26)+nrow(migratory_proportion_df_week_1), ncol = 2))
colnames(nspecies_hist) <- c("value", "type")
nspecies_hist[c(1:nrow(migratory_proportion_df_week_1)), 1] <- migratory_proportion_df_week_1$nspecies
nspecies_hist[c((nrow(migratory_proportion_df_week_1) + 1)):(nrow(nspecies_hist)),1] <- migratory_proportion_df_week_26$nspecies
nspecies_hist[c(1:nrow(migratory_proportion_df_week_1)), 2] <- "Non-Breeding"
nspecies_hist[c((nrow(migratory_proportion_df_week_1) + 1)):(nrow(nspecies_hist)),2] <- "Breeding"
tail(nspecies_hist)

ggplot(nspecies_hist, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.5) +
  xlab("Species Richness") +
  ylab("Density") + 
  guides(fill=guide_legend(title="Season")) +
  theme(legend.position = "bottom", text = element_text(size = 20))


### add migratory proportion to the other data frames

for(i in 1:52){
  assign(paste0("migratory_proportion_df_week", i), read.csv(paste0("output/community_data/migratory_proportion_week_", i, ".csv")))
} #takes about a half hour

tail(all_combined_week_NRIs_df)
head(migratory_proportion_df_week_1)
migratory_proportion_df_week_1[migratory_proportion_df_week_1$community == all_combined_week_NRIs_df[1,"X"],"proportion"]
migratory_proportion_df_week2[migratory_proportion_df_week2$community == all_combined_week_NRIs_df[1,"X"],"proportion"]

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

all_combined_week_NRIs_df[1,"X"]
all_combined_week_NRIs_df$proportion <- rep(NA, nrow(all_combined_week_NRIs_df))
# for(i in 1:nrow(all_combined_week_NRIs_df)){ #changed value based on where it left off
#   all_combined_week_NRIs_df[i,"proportion"] <- all_mig_props[[all_combined_week_NRIs_df[i,".id"]]][all_mig_props[[all_combined_week_NRIs_df[i,".id"]]]$community == all_combined_week_NRIs_df[i,"X"],"proportion"]
# }
#  warning: takes a long time (~3 days)
# write.csv(all_combined_week_NRIs_df, "output/tables/combined_NRI_df.csv")
all_combined_week_NRIs_df <- read.csv("output/tables/combined_NRI_df.csv")

#replicate with NTI
all_combined_week_NTIs_df$proportion <- rep(NA, nrow(all_combined_week_NTIs_df))
# for(i in 3606:nrow(all_combined_week_NTIs_df)){ #changed value based on where it left off
#   all_combined_week_NTIs_df[i,"proportion"] <- all_mig_props[[all_combined_week_NTIs_df[i,".id"]]][all_mig_props[[all_combined_week_NTIs_df[i,".id"]]]$community == all_combined_week_NTIs_df[i,"X"],"proportion"]
# }
# write.csv(all_combined_week_NTIs_df, "output/tables/combined_NTI_df.csv")
all_combined_week_NTIs_df <- read.csv("output/tables/combined_NTI_df.csv")
# warning: takes a long time (~3 days)




