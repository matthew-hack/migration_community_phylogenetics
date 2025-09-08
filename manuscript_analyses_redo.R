# 5. Manuscript data analysis
# Generates primary data analyses included in the manuscript. 
# Sets up seasonal comparisons to model seasonal turnover in phylogenetic community metrics against other factors.
# Generates figures and tables relating to seasonal turnover in phylogenetic community metrics.

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

# load the mapping parameters to generate necessary CRS
path <- get_species_path("example_data")
fac_parameters <- load_fac_map_parameters(path)
crs <- fac_parameters$custom_projection

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

##### CREATE A HISTOGRAM WITH DISTRIBUTIONS OF PHYLOGENETIC DISTANCE MEANS ACROSS ALL BCRS (SUPPLEMENTAL FIGURE 1) #####

#load in ecoregional community data for week 26 (summer)
for(i in 5:37){
  assign(paste0("bcr_", i, "_week_26"), read.csv(paste0("output/bcr_community_data/26/community_26_bcr_", i, ".csv")))
} #takes about a half hour

#load in ecoregional community data for week 1 (winter)
for(i in 5:37){
  assign(paste0("bcr_", i, "_week_1"), read.csv(paste0("output/bcr_community_data/1/community_1_bcr_", i, ".csv")))
} #takes about a half hour

### process summer data
#create list of all week 26 community matrices
bcr_week_26_list <- list(bcr_5_week_26, bcr_6_week_26, bcr_7_week_26, bcr_8_week_26, bcr_9_week_26, bcr_10_week_26,
                         bcr_11_week_26, bcr_12_week_26, bcr_13_week_26, bcr_14_week_26, bcr_15_week_26, bcr_16_week_26,
                         bcr_17_week_26, bcr_18_week_26, bcr_19_week_26, bcr_20_week_26, bcr_21_week_26, bcr_22_week_26,
                         bcr_23_week_26, bcr_24_week_26, bcr_25_week_26, bcr_26_week_26, bcr_27_week_26, bcr_28_week_26,
                         bcr_29_week_26, bcr_30_week_26, bcr_31_week_26, bcr_32_week_26, bcr_33_week_26, bcr_34_week_26,
                         bcr_35_week_26, bcr_36_week_26, bcr_37_week_26)

#generate c scores for all species combinations (denotes whether species combinations are found in sympatry more or less often than random chance)
for(i in 5:37){
  df <- as.data.frame(bcr_week_26_list[i-4])
  df <- df[sample(nrow(df), nrow(df)/100, replace=FALSE),]
  ecospot <- ecospat.Cscore(df, 2, outpath="output/tables")
  assign(paste0("c_scores_bcr_", i, "_week_26"), read.delim("output/tables/Cscores.txt"))
}

#create list of all c score matrices
c_scores_bcr_week_26_list <- list(c_scores_bcr_5_week_26, c_scores_bcr_6_week_26, c_scores_bcr_7_week_26, c_scores_bcr_8_week_26, c_scores_bcr_9_week_26, c_scores_bcr_10_week_26,
                                  c_scores_bcr_11_week_26, c_scores_bcr_12_week_26, c_scores_bcr_13_week_26, c_scores_bcr_14_week_26, c_scores_bcr_15_week_26, c_scores_bcr_16_week_26,
                                  c_scores_bcr_17_week_26, c_scores_bcr_18_week_26, c_scores_bcr_19_week_26, c_scores_bcr_20_week_26, c_scores_bcr_21_week_26, c_scores_bcr_22_week_26,
                                  c_scores_bcr_23_week_26, c_scores_bcr_24_week_26, c_scores_bcr_25_week_26, c_scores_bcr_26_week_26, c_scores_bcr_27_week_26, c_scores_bcr_28_week_26,
                                  c_scores_bcr_29_week_26, c_scores_bcr_30_week_26, c_scores_bcr_31_week_26, c_scores_bcr_32_week_26, c_scores_bcr_33_week_26, c_scores_bcr_34_week_26,
                                  c_scores_bcr_35_week_26, c_scores_bcr_36_week_26, c_scores_bcr_37_week_26)

#generate cleaned phylogenetic trees and community matrices for each BCR
for(i in 5:37){
  assign(paste0("passerine.picCleanTree.bcr", i, "week26"), match.phylo.comm(phy = passerine.tree,
                                                                             comm = as.data.frame(bcr_week_26_list[i-4]))$phy)
  assign(paste0("passerine.picCleanComm.bcr", i, "week26"), match.phylo.comm(phy = passerine.tree,
                                                                             comm = as.data.frame(bcr_week_26_list[i-4]))$comm)
}

#make a list of cleaned trees
passerine.picCleanTree.bcrweek_26_list <- list(passerine.picCleanTree.bcr5week26, passerine.picCleanTree.bcr6week26, passerine.picCleanTree.bcr7week26, passerine.picCleanTree.bcr8week26, passerine.picCleanTree.bcr9week26, passerine.picCleanTree.bcr10week26,
                                               passerine.picCleanTree.bcr11week26, passerine.picCleanTree.bcr12week26, passerine.picCleanTree.bcr13week26, passerine.picCleanTree.bcr14week26, passerine.picCleanTree.bcr15week26, passerine.picCleanTree.bcr16week26,
                                               passerine.picCleanTree.bcr17week26, passerine.picCleanTree.bcr18week26, passerine.picCleanTree.bcr19week26, passerine.picCleanTree.bcr20week26, passerine.picCleanTree.bcr21week26, passerine.picCleanTree.bcr22week26,
                                               passerine.picCleanTree.bcr23week26, passerine.picCleanTree.bcr24week26, passerine.picCleanTree.bcr25week26, passerine.picCleanTree.bcr26week26, passerine.picCleanTree.bcr27week26, passerine.picCleanTree.bcr28week26,
                                               passerine.picCleanTree.bcr29week26, passerine.picCleanTree.bcr30week26, passerine.picCleanTree.bcr31week26, passerine.picCleanTree.bcr32week26, passerine.picCleanTree.bcr33week26, passerine.picCleanTree.bcr34week26,
                                               passerine.picCleanTree.bcr35week26, passerine.picCleanTree.bcr36week26, passerine.picCleanTree.bcr37week26)

# create phylogenetic distance matrices and make a list of those matrices
for(i in 5:37){
  assign(paste0("passerine.cophenDist.bcr", i, "week26"), cophenetic.phylo(passerine.picCleanTree.bcrweek_26_list[[i-4]]))
}
passerine.cophenDist.week26_list <- list(passerine.cophenDist.bcr5week26, passerine.cophenDist.bcr6week26, passerine.cophenDist.bcr7week26, passerine.cophenDist.bcr8week26, passerine.cophenDist.bcr9week26, passerine.cophenDist.bcr10week26,
                                         passerine.cophenDist.bcr11week26, passerine.cophenDist.bcr12week26, passerine.cophenDist.bcr13week26, passerine.cophenDist.bcr14week26, passerine.cophenDist.bcr15week26, passerine.cophenDist.bcr16week26,
                                         passerine.cophenDist.bcr17week26, passerine.cophenDist.bcr18week26, passerine.cophenDist.bcr19week26, passerine.cophenDist.bcr20week26, passerine.cophenDist.bcr21week26, passerine.cophenDist.bcr22week26,
                                         passerine.cophenDist.bcr23week26, passerine.cophenDist.bcr24week26, passerine.cophenDist.bcr25week26, passerine.cophenDist.bcr26week26, passerine.cophenDist.bcr27week26, passerine.cophenDist.bcr28week26,
                                         passerine.cophenDist.bcr29week26, passerine.cophenDist.bcr30week26, passerine.cophenDist.bcr31week26, passerine.cophenDist.bcr32week26, passerine.cophenDist.bcr33week26, passerine.cophenDist.bcr34week26,
                                         passerine.cophenDist.bcr35week26, passerine.cophenDist.bcr36week26, passerine.cophenDist.bcr37week26)

#populate phylogenetic distances for each combination of species in each BCR
for(i in 5:37){
  c_scores_bcr_week_26_list[[i-4]]$phylogenetic_distance <- rep(NA, nrow(c_scores_bcr_week_26_list[[i-4]]))
  for(j in 1:nrow(c_scores_bcr_week_26_list[[i-4]])){
    c_scores_bcr_week_26_list[[i-4]][j, "phylogenetic_distance"] <- passerine.cophenDist.week26_list[[i-4]][c_scores_bcr_week_26_list[[i-4]][j,"Sp1"],
                                                                                                            c_scores_bcr_week_26_list[[i-4]][j,"Sp2"]]
  }
}

#make a vector of mean NRI value within each BCR
week_26_mean_NRI_values <- rep(NA, length(c_scores_bcr_week_26_list))
for(i in 1:length(c_scores_bcr_week_26_list)){
  week_26_mean_NRI_values[i] <- mean(c_scores_bcr_week_26_list[[i]]$phylogenetic_distance)
}

#calculate nearest neighbors for each matrix
#generate species lists for each BCR in week 26
for(i in 5:37){
  assign(paste0("bcr", i, "week26_species"), colnames(bcr_week_26_list[[i-4]]))
}

#generate list of species lists
week26_species_lists <- list(bcr5week26_species, bcr6week26_species, bcr7week26_species, bcr8week26_species, bcr9week26_species, bcr10week26_species,
                             bcr11week26_species, bcr12week26_species, bcr13week26_species, bcr14week26_species, bcr15week26_species, bcr16week26_species,
                             bcr17week26_species, bcr18week26_species, bcr19week26_species, bcr20week26_species, bcr21week26_species, bcr22week26_species,
                             bcr23week26_species, bcr24week26_species, bcr25week26_species, bcr26week26_species, bcr27week26_species, bcr28week26_species,
                             bcr29week26_species, bcr30week26_species, bcr31week26_species, bcr32week26_species, bcr33week26_species, bcr34week26_species,
                             bcr35week26_species, bcr36week26_species, bcr37week26_species)

#generate phylogenetic distance matrices
for(i in 5:37){
  assign(paste0("bcr", i, "week26_min_phy_distances"), rep(NA, length(week26_species_lists[[i-4]])))
}

#generate list of phylogenetic distance matrices
week26_min_phy_distances_lists <- list(bcr5week26_min_phy_distances, bcr6week26_min_phy_distances, bcr7week26_min_phy_distances, bcr8week26_min_phy_distances, bcr9week26_min_phy_distances, bcr10week26_min_phy_distances,
                                       bcr11week26_min_phy_distances, bcr12week26_min_phy_distances, bcr13week26_min_phy_distances, bcr14week26_min_phy_distances, bcr15week26_min_phy_distances, bcr16week26_min_phy_distances,
                                       bcr17week26_min_phy_distances, bcr18week26_min_phy_distances, bcr19week26_min_phy_distances, bcr20week26_min_phy_distances, bcr21week26_min_phy_distances, bcr22week26_min_phy_distances,
                                       bcr23week26_min_phy_distances, bcr24week26_min_phy_distances, bcr25week26_min_phy_distances, bcr26week26_min_phy_distances, bcr27week26_min_phy_distances, bcr28week26_min_phy_distances,
                                       bcr29week26_min_phy_distances, bcr30week26_min_phy_distances, bcr31week26_min_phy_distances, bcr32week26_min_phy_distances, bcr33week26_min_phy_distances, bcr34week26_min_phy_distances,
                                       bcr35week26_min_phy_distances, bcr36week26_min_phy_distances, bcr37week26_min_phy_distances)

#populate matrices of minimum phylogenetic distances
for(i in 5:37){
  for(j in 1:length(week26_min_phy_distances_lists[[i-4]])){
    week26_min_phy_distances_lists[[i-4]][j] <- min(c_scores_bcr_week_26_list[[i-4]][c_scores_bcr_week_26_list[[i-4]]$Sp1 == week26_species_lists[[i-4]][j] |
                                                                                       c_scores_bcr_week_26_list[[i-4]]$Sp2 == week26_species_lists[[i-4]][j],
                                                                                     "phylogenetic_distance"])
  }
}

#generate a vector of mean NTI value within each BCR
week_26_mean_NTI_values <- rep(NA, length(c_scores_bcr_week_26_list))
for(i in 1:length(c_scores_bcr_week_26_list)){
  week_26_mean_NTI_values[i] <- mean(week26_min_phy_distances_lists[[i]])
}

### replicate process for winter data

#create list of all week 1 community matrices
bcr_week_1_list <- list(bcr_5_week_1, bcr_6_week_1, bcr_7_week_1, bcr_8_week_1, bcr_9_week_1, bcr_10_week_1,
                        bcr_11_week_1, bcr_12_week_1, bcr_13_week_1, bcr_14_week_1, bcr_15_week_1, bcr_16_week_1,
                        bcr_17_week_1, bcr_18_week_1, bcr_19_week_1, bcr_20_week_1, bcr_21_week_1, bcr_22_week_1,
                        bcr_23_week_1, bcr_24_week_1, bcr_25_week_1, bcr_26_week_1, bcr_27_week_1, bcr_28_week_1,
                        bcr_29_week_1, bcr_30_week_1, bcr_31_week_1, bcr_32_week_1, bcr_33_week_1, bcr_34_week_1,
                        bcr_35_week_1, bcr_36_week_1, bcr_37_week_1)

#generate c scores for all species combinations (denotes whether species combinations are found in sympatry more or less often than random chance)
for(i in 5:37){
  df <- as.data.frame(bcr_week_1_list[i-4])
  df <- df[sample(nrow(df), nrow(df)/100, replace=FALSE),]
  ecospot <- ecospat.Cscore(df, 2, outpath="output/tables")
  assign(paste0("c_scores_bcr_", i, "_week_1"), read.delim("output/tables/Cscores.txt"))
}

#create list of all c score matrices
c_scores_bcr_week_1_list <- list(c_scores_bcr_5_week_1, c_scores_bcr_6_week_1, c_scores_bcr_7_week_1, c_scores_bcr_8_week_1, c_scores_bcr_9_week_1, c_scores_bcr_10_week_1,
                                 c_scores_bcr_11_week_1, c_scores_bcr_12_week_1, c_scores_bcr_13_week_1, c_scores_bcr_14_week_1, c_scores_bcr_15_week_1, c_scores_bcr_16_week_1,
                                 c_scores_bcr_17_week_1, c_scores_bcr_18_week_1, c_scores_bcr_19_week_1, c_scores_bcr_20_week_1, c_scores_bcr_21_week_1, c_scores_bcr_22_week_1,
                                 c_scores_bcr_23_week_1, c_scores_bcr_24_week_1, c_scores_bcr_25_week_1, c_scores_bcr_26_week_1, c_scores_bcr_27_week_1, c_scores_bcr_28_week_1,
                                 c_scores_bcr_29_week_1, c_scores_bcr_30_week_1, c_scores_bcr_31_week_1, c_scores_bcr_32_week_1, c_scores_bcr_33_week_1, c_scores_bcr_34_week_1,
                                 c_scores_bcr_35_week_1, c_scores_bcr_36_week_1, c_scores_bcr_37_week_1)

#generate cleaned phylogenetic trees and community matrices for each BCR
for(i in 5:37){
  assign(paste0("passerine.picCleanTree.bcr", i, "week1"), match.phylo.comm(phy = passerine.tree,
                                                                            comm = as.data.frame(bcr_week_1_list[i-4]))$phy)
  assign(paste0("passerine.picCleanComm.bcr", i, "week1"), match.phylo.comm(phy = passerine.tree,
                                                                            comm = as.data.frame(bcr_week_1_list[i-4]))$comm)
}

#make a list of cleaned trees
passerine.picCleanTree.bcrweek_1_list <- list(passerine.picCleanTree.bcr5week1, passerine.picCleanTree.bcr6week1, passerine.picCleanTree.bcr7week1, passerine.picCleanTree.bcr8week1, passerine.picCleanTree.bcr9week1, passerine.picCleanTree.bcr10week1,
                                              passerine.picCleanTree.bcr11week1, passerine.picCleanTree.bcr12week1, passerine.picCleanTree.bcr13week1, passerine.picCleanTree.bcr14week1, passerine.picCleanTree.bcr15week1, passerine.picCleanTree.bcr16week1,
                                              passerine.picCleanTree.bcr17week1, passerine.picCleanTree.bcr18week1, passerine.picCleanTree.bcr19week1, passerine.picCleanTree.bcr20week1, passerine.picCleanTree.bcr21week1, passerine.picCleanTree.bcr22week1,
                                              passerine.picCleanTree.bcr23week1, passerine.picCleanTree.bcr24week1, passerine.picCleanTree.bcr25week1, passerine.picCleanTree.bcr26week1, passerine.picCleanTree.bcr27week1, passerine.picCleanTree.bcr28week1,
                                              passerine.picCleanTree.bcr29week1, passerine.picCleanTree.bcr30week1, passerine.picCleanTree.bcr31week1, passerine.picCleanTree.bcr32week1, passerine.picCleanTree.bcr33week1, passerine.picCleanTree.bcr34week1,
                                              passerine.picCleanTree.bcr35week1, passerine.picCleanTree.bcr36week1, passerine.picCleanTree.bcr37week1)

# create phylogenetic distance matrices and make a list of those matrices
for(i in 5:37){
  assign(paste0("passerine.cophenDist.bcr", i, "week1"), cophenetic.phylo(passerine.picCleanTree.bcrweek_1_list[[i-4]]))
}

#make a list of cleaned trees
passerine.cophenDist.week1_list <- list(passerine.cophenDist.bcr5week1, passerine.cophenDist.bcr6week1, passerine.cophenDist.bcr7week1, passerine.cophenDist.bcr8week1, passerine.cophenDist.bcr9week1, passerine.cophenDist.bcr10week1,
                                        passerine.cophenDist.bcr11week1, passerine.cophenDist.bcr12week1, passerine.cophenDist.bcr13week1, passerine.cophenDist.bcr14week1, passerine.cophenDist.bcr15week1, passerine.cophenDist.bcr16week1,
                                        passerine.cophenDist.bcr17week1, passerine.cophenDist.bcr18week1, passerine.cophenDist.bcr19week1, passerine.cophenDist.bcr20week1, passerine.cophenDist.bcr21week1, passerine.cophenDist.bcr22week1,
                                        passerine.cophenDist.bcr23week1, passerine.cophenDist.bcr24week1, passerine.cophenDist.bcr25week1, passerine.cophenDist.bcr26week1, passerine.cophenDist.bcr27week1, passerine.cophenDist.bcr28week1,
                                        passerine.cophenDist.bcr29week1, passerine.cophenDist.bcr30week1, passerine.cophenDist.bcr31week1, passerine.cophenDist.bcr32week1, passerine.cophenDist.bcr33week1, passerine.cophenDist.bcr34week1,
                                        passerine.cophenDist.bcr35week1, passerine.cophenDist.bcr36week1, passerine.cophenDist.bcr37week1)

#populate phylogenetic distances for each combination of species in each BCR
for(i in 5:37){
  c_scores_bcr_week_1_list[[i-4]]$phylogenetic_distance <- rep(NA, nrow(c_scores_bcr_week_1_list[[i-4]]))
  for(j in 1:nrow(c_scores_bcr_week_1_list[[i-4]])){
    c_scores_bcr_week_1_list[[i-4]][j, "phylogenetic_distance"] <- passerine.cophenDist.week1_list[[i-4]][c_scores_bcr_week_1_list[[i-4]][j,"Sp1"],
                                                                                                          c_scores_bcr_week_1_list[[i-4]][j,"Sp2"]]
  }
}

#make a vector of mean NRI value within each BCR
week_1_mean_NRI_values <- rep(NA, length(c_scores_bcr_week_1_list))
for(i in 1:length(c_scores_bcr_week_1_list)){
  week_1_mean_NRI_values[i] <- mean(c_scores_bcr_week_1_list[[i]]$phylogenetic_distance)
}

#calculate nearest neighbors for each matrix
#generate species lists for each BCR in week 1
for(i in 5:37){
  assign(paste0("bcr", i, "week1_species"), colnames(bcr_week_1_list[[i-4]]))
}

#generate list of species lists
week1_species_lists <- list(bcr5week1_species, bcr6week1_species, bcr7week1_species, bcr8week1_species, bcr9week1_species, bcr10week1_species,
                            bcr11week1_species, bcr12week1_species, bcr13week1_species, bcr14week1_species, bcr15week1_species, bcr16week1_species,
                            bcr17week1_species, bcr18week1_species, bcr19week1_species, bcr20week1_species, bcr21week1_species, bcr22week1_species,
                            bcr23week1_species, bcr24week1_species, bcr25week1_species, bcr26week1_species, bcr27week1_species, bcr28week1_species,
                            bcr29week1_species, bcr30week1_species, bcr31week1_species, bcr32week1_species, bcr33week1_species, bcr34week1_species,
                            bcr35week1_species, bcr36week1_species, bcr37week1_species)

#generate phylogenetic distance matrices
for(i in 5:37){
  assign(paste0("bcr", i, "week1_min_phy_distances"), rep(NA, length(week1_species_lists[[i-4]])))
}

#generate list of phylogenetic distance matrices
week1_min_phy_distances_lists <- list(bcr5week1_min_phy_distances, bcr6week1_min_phy_distances, bcr7week1_min_phy_distances, bcr8week1_min_phy_distances, bcr9week1_min_phy_distances, bcr10week1_min_phy_distances,
                                      bcr11week1_min_phy_distances, bcr12week1_min_phy_distances, bcr13week1_min_phy_distances, bcr14week1_min_phy_distances, bcr15week1_min_phy_distances, bcr16week1_min_phy_distances,
                                      bcr17week1_min_phy_distances, bcr18week1_min_phy_distances, bcr19week1_min_phy_distances, bcr20week1_min_phy_distances, bcr21week1_min_phy_distances, bcr22week1_min_phy_distances,
                                      bcr23week1_min_phy_distances, bcr24week1_min_phy_distances, bcr25week1_min_phy_distances, bcr26week1_min_phy_distances, bcr27week1_min_phy_distances, bcr28week1_min_phy_distances,
                                      bcr29week1_min_phy_distances, bcr30week1_min_phy_distances, bcr31week1_min_phy_distances, bcr32week1_min_phy_distances, bcr33week1_min_phy_distances, bcr34week1_min_phy_distances,
                                      bcr35week1_min_phy_distances, bcr36week1_min_phy_distances, bcr37week1_min_phy_distances)

#populate matrices of minimum phylogenetic distances
for(i in 5:37){
  for(j in 1:length(week1_min_phy_distances_lists[[i-4]])){
    week1_min_phy_distances_lists[[i-4]][j] <- min(c_scores_bcr_week_1_list[[i-4]][c_scores_bcr_week_1_list[[i-4]]$Sp1 == week1_species_lists[[i-4]][j] |
                                                                                     c_scores_bcr_week_1_list[[i-4]]$Sp2 == week1_species_lists[[i-4]][j],
                                                                                   "phylogenetic_distance"])
  }
}

#generate a vector of mean NTI value within each BCR
week_1_mean_NTI_values <- rep(NA, length(c_scores_bcr_week_1_list))
for(i in 1:length(c_scores_bcr_week_1_list)){
  week_1_mean_NTI_values[i] <- mean(week1_min_phy_distances_lists[[i]])
}

### create density plots for distributions of mean NRI and mean NTI for weeks 1 and 26 (SF1)
#process data for week 26
#generate data frame for plotting mean NRI and mean NTI data for week 26
mean_NRI_NTI_hist_26 <- data.frame(matrix(NA, nrow = 2*length(week_26_mean_NRI_values), ncol = 2))
colnames(mean_NRI_NTI_hist_26) <- c("value", "type")
#populate mean NRI values
mean_NRI_NTI_hist_26[c(1:length(week_26_mean_NRI_values)), "value"] <- week_26_mean_NRI_values
#populate mean NTI values
mean_NRI_NTI_hist_26[c((length(week_26_mean_NRI_values) + 1)):(nrow(mean_NRI_NTI_hist_26)), "value"] <- week_26_mean_NTI_values
#add labels
mean_NRI_NTI_hist_26[c(1:length(week_26_mean_NRI_values)), "type"] <- "NRI"
mean_NRI_NTI_hist_26[c((length(week_26_mean_NRI_values) + 1)):(nrow(mean_NRI_NTI_hist_26)), "type"] <- "NTI"

#create gg object depicting week 26 mean NRI and NTI distributions 
week26density <- ggplot(mean_NRI_NTI_hist_26, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("Phylogenetic Distance") +
  ylab("Density") +
  xlim(0,100) +
  guides(fill=guide_legend(title="Metric")) +
  ggtitle("Week 26") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25))

#process data for week 1
#generate data frame for plotting mean NRI and mean NTI data for week 1
mean_NRI_NTI_hist_1 <- data.frame(matrix(NA, nrow = 2*length(week_1_mean_NRI_values), ncol = 2))
colnames(mean_NRI_NTI_hist_1) <- c("value", "type")
#populate mean NRI values
mean_NRI_NTI_hist_1[c(1:length(week_1_mean_NRI_values)), "value"] <- week_1_mean_NRI_values
#populate mean NTI values
mean_NRI_NTI_hist_1[c((length(week_1_mean_NRI_values) + 1)):(nrow(mean_NRI_NTI_hist_1)), "value"] <- week_1_mean_NTI_values
#add labels
mean_NRI_NTI_hist_1[c(1:length(week_1_mean_NRI_values)), "type"] <- "NRI"
mean_NRI_NTI_hist_1[c((length(week_1_mean_NRI_values) + 1)):(nrow(mean_NRI_NTI_hist_1)), "type"] <- "NTI"

#create gg object depicting week 1 mean NRI and NTI distributions 
week1density <- ggplot(mean_NRI_NTI_hist_1, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("Phylogenetic Distance") +
  ylab("Density") +
  xlim(0,100) +
  guides(fill=guide_legend(title="Metric")) +
  ggtitle("Week 1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25))

#merge dataframes for the two weeks and label
both_weeks_density <- dplyr::bind_rows(mean_NRI_NTI_hist_26, mean_NRI_NTI_hist_1, .id = "week")
both_weeks_density[both_weeks_density$week == "2","week"] <- "Week 1"
both_weeks_density[both_weeks_density$week == "1","week"] <- "Week 26"

#arrange two weeks' density plots into one plot
ggarrange(week1density, week26density, ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")

#generate Supplemental Figure 1
pdf(file = 'output/figures/mean_NRI_NTI_values_density.pdf')
ggplot(both_weeks_density, aes(x=value, group=type, fill=type)) +
  geom_density(alpha=0.67) +
  xlab("Phylogenetic Distance") +
  ylab("Density") +
  xlim(0,100) +
  guides(fill=guide_legend(title="Metric")) +
  ggtitle("Mean NRI vs. NTI Phylogenetic Distances") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~week) +
  scale_fill_manual(values=c("red","blue")) +
  theme(legend.position = "bottom", text = element_text(size = 25))
dev.off()

##### PROCESS WEEK 1 AND WEEK 26 DATASETS FOR SEASONAL TURNOVER ANALYSES #####

### GO BACK: change redo names
#reread in NRI and NTI ecoregional data generated from weeks 1 and weeks 26
combined_week_26_NRI <- read.csv("output/bcr_tables/26/combined_NRI_26_redo.csv")
combined_week_26_NTI <- read.csv("output/bcr_tables/26/combined_NTI_26_redo.csv")
combined_week_1_NRI <- read.csv("output/bcr_tables/1/combined_NRI_1_redo.csv")
combined_week_1_NTI <- read.csv("output/bcr_tables/1/combined_NTI_1_redo.csv")

#rename community column as "community" if it has been renamed as "X" (skip if not necessary)
combined_week_1_NTI <- combined_week_1_NTI %>%
  rename(community = X)
combined_week_1_NRI <- combined_week_1_NRI %>%
  rename(community = X)
combined_week_26_NTI <- combined_week_26_NTI %>%
  rename(community = X)
combined_week_26_NRI <- combined_week_26_NRI %>%
  rename(community = X)

### GO BACK: this should have been done in past steps and so therefore isn't necessary
#load migratory proportion files for weeks 1 and 26
migratory_proportion_df_week_1 <- read.csv("output/community_data/migratory_proportion_week_1.csv")
migratory_proportion_df_week_26 <- read.csv("output/community_data/migratory_proportion_week_26.csv")

#add migratory proportion field to each dataframe
combined_week_1_NRI$proportion <- rep(NA, nrow(combined_week_1_NRI))
combined_week_1_NTI$proportion <- rep(NA, nrow(combined_week_1_NTI))
combined_week_26_NRI$proportion <- rep(NA, nrow(combined_week_26_NRI))
combined_week_26_NTI$proportion <- rep(NA, nrow(combined_week_26_NTI))

#populate proportion fields in each dataframe
for(i in 1:nrow(combined_week_1_NRI)){
  combined_week_1_NRI[i,"proportion"] <- migratory_proportion_df_week_1[migratory_proportion_df_week_1$community == combined_week_1_NRI[i,"community"],"proportion"]
  combined_week_1_NTI[i,"proportion"] <- migratory_proportion_df_week_1[migratory_proportion_df_week_1$community == combined_week_1_NTI[i,"community"],"proportion"]
}
for(i in 1:nrow(combined_week_26_NRI)){
  combined_week_26_NRI[i,"proportion"] <- migratory_proportion_df_week_26[migratory_proportion_df_week_26$community == combined_week_26_NRI[i,"community"],"proportion"]
  combined_week_26_NTI[i,"proportion"] <- migratory_proportion_df_week_26[migratory_proportion_df_week_26$community == combined_week_26_NTI[i,"community"],"proportion"]
}

#convert x and y values to latitude and longitude
combined_week_1_NRI$latitude <- sf_project(from = crs, to = st_crs(4326), combined_week_1_NRI[,c("x","y")])[,2]
combined_week_1_NRI$longitude <- sf_project(from = crs, to = st_crs(4326), combined_week_1_NRI[,c("x","y")])[,1]
combined_week_1_NTI$latitude <- sf_project(from = crs, to = st_crs(4326), combined_week_1_NTI[,c("x","y")])[,2]
combined_week_1_NTI$longitude <- sf_project(from = crs, to = st_crs(4326), combined_week_1_NTI[,c("x","y")])[,1]
combined_week_26_NRI$latitude <- sf_project(from = crs, to = st_crs(4326), combined_week_26_NRI[,c("x","y")])[,2]
combined_week_26_NRI$longitude <- sf_project(from = crs, to = st_crs(4326), combined_week_26_NRI[,c("x","y")])[,1]
combined_week_26_NTI$latitude <- sf_project(from = crs, to = st_crs(4326), combined_week_26_NTI[,c("x","y")])[,2]
combined_week_26_NTI$longitude <- sf_project(from = crs, to = st_crs(4326), combined_week_26_NTI[,c("x","y")])[,1]

##### CALCULATE SEASONAL TURNOVER VALUES IN PHYLOGENETIC COMMUNITY METRICS FOR ECOREGIONAL COMMUNITIES #####

### PROCESS NTI DATA

#create data frames representing weeks 1 and 26 to initalize seasonal turnover dataframe
NTI_week_1 <- rep(NA, nrow(combined_week_1_NTI))
NTI_week_26 <- rep(NA, nrow(combined_week_26_NTI))
#merge into seasonal turnover into one dataframe
NTI_turnover_ratios <- data.frame(NTI_week_1, NTI_week_26)
#add fields for species richness
NTI_turnover_ratios$richness_week_1 <- rep(NA, nrow(combined_week_1_NTI))
NTI_turnover_ratios$richness_week_26 <- rep(NA, nrow(combined_week_1_NTI))
#add field for seasonal turnover in NTI
NTI_turnover_ratios$NTI_turnover <- rep(NA, nrow(combined_week_1_NTI))
#add field for seasonal turnover in species richness
NTI_turnover_ratios$richness_turnover <- rep(NA, nrow(combined_week_1_NTI))
#add fields for BCR and geographic coordinates
NTI_turnover_ratios$BCR <- rep(NA, nrow(combined_week_1_NTI))
NTI_turnover_ratios$latitude <- rep(NA, nrow(combined_week_1_NTI))
NTI_turnover_ratios$longitude <- rep(NA, nrow(combined_week_1_NTI))
#add fields for migratory proportion
NTI_turnover_ratios$mig_prop_week_1 <- rep(NA, nrow(combined_week_1_NTI))
NTI_turnover_ratios$mig_prop_week_26 <- rep(NA, nrow(combined_week_1_NTI))
#add field for seasonal turnover in migratory proportion
NTI_turnover_ratios$mig_prop_turnover <- rep(NA, nrow(combined_week_1_NTI))
#add field for individual community identifiers
NTI_turnover_ratios$community <- rep(NA, nrow(combined_week_1_NTI))

#populate dataframe generated above

#add NTI values
NTI_turnover_ratios$NTI_week_1 <- combined_week_1_NTI$NTI
NTI_turnover_ratios$NTI_week_26 <- combined_week_26_NTI$NTI
#add species richness values
NTI_turnover_ratios$richness_week_1 <- combined_week_1_NTI$ntaxa
NTI_turnover_ratios$richness_week_26 <- combined_week_26_NTI$ntaxa
#calculate turnover in NTI and species richness
NTI_turnover_ratios$NTI_turnover <- NTI_turnover_ratios$NTI_week_26 - NTI_turnover_ratios$NTI_week_1
NTI_turnover_ratios$richness_turnover <- NTI_turnover_ratios$richness_week_26 - NTI_turnover_ratios$richness_week_1
#populate BCR and geographic coordinates
NTI_turnover_ratios$BCR <- combined_week_1_NTI$BCR
NTI_turnover_ratios$latitude <- combined_week_1_NTI$latitude
NTI_turnover_ratios$longitude <- combined_week_1_NTI$longitude
#add migratory proportion values
NTI_turnover_ratios$mig_prop_week_1 <- combined_week_1_NTI$proportion
NTI_turnover_ratios$mig_prop_week_26 <- combined_week_26_NTI$proportion
#calculate turnover in migratory proportion
NTI_turnover_ratios$mig_prop_turnover <- NTI_turnover_ratios$mig_prop_week_26 - NTI_turnover_ratios$mig_prop_week_1
#add individual community identifier
NTI_turnover_ratios$community <- combined_week_1_NTI$community

### PROCESS NRI DATA

#create data frames representing weeks 1 and 26 to initalize seasonal turnover dataframe
NRI_week_1 <- rep(NA, nrow(combined_week_1_NRI))
NRI_week_26 <- rep(NA, nrow(combined_week_26_NRI))
#merge into seasonal turnover into one dataframe
NRI_turnover_ratios <- data.frame(NRI_week_1, NRI_week_26)
#add fields for species richness
NRI_turnover_ratios$richness_week_1 <- rep(NA, nrow(combined_week_1_NRI))
NRI_turnover_ratios$richness_week_26 <- rep(NA, nrow(combined_week_1_NRI))
#add field for seasonal turnover in NRI
NRI_turnover_ratios$NRI_turnover <- rep(NA, nrow(combined_week_1_NRI))
#add field for seasonal turnover in species richness
NRI_turnover_ratios$richness_turnover <- rep(NA, nrow(combined_week_1_NRI))
#add fields for BCR and geographic coordinates
NRI_turnover_ratios$BCR <- rep(NA, nrow(combined_week_1_NRI))
NRI_turnover_ratios$latitude <- rep(NA, nrow(combined_week_1_NRI))
NRI_turnover_ratios$longitude <- rep(NA, nrow(combined_week_1_NRI))
#add fields for migratory proportion
NRI_turnover_ratios$mig_prop_week_1 <- rep(NA, nrow(combined_week_1_NRI))
NRI_turnover_ratios$mig_prop_week_26 <- rep(NA, nrow(combined_week_1_NRI))
#add field for seasonal turnover in migratory proportion
NRI_turnover_ratios$mig_prop_turnover <- rep(NA, nrow(combined_week_1_NRI))
#add field for individual community identifiers
NRI_turnover_ratios$community <- rep(NA, nrow(combined_week_1_NRI))

#populate dataframe generated above

#add NRI values
NRI_turnover_ratios$NRI_week_1 <- combined_week_1_NRI$NRI
NRI_turnover_ratios$NRI_week_26 <- combined_week_26_NRI$NRI
#add species richness values
NRI_turnover_ratios$richness_week_1 <- combined_week_1_NRI$ntaxa
NRI_turnover_ratios$richness_week_26 <- combined_week_26_NRI$ntaxa
#calculate turnover in NRI and species richness
NRI_turnover_ratios$NRI_turnover <- NRI_turnover_ratios$NRI_week_26 - NRI_turnover_ratios$NRI_week_1
NRI_turnover_ratios$richness_turnover <- NRI_turnover_ratios$richness_week_26 - NRI_turnover_ratios$richness_week_1
#populate BCR and geographic coordinates
NRI_turnover_ratios$BCR <- combined_week_1_NRI$BCR
NRI_turnover_ratios$latitude <- combined_week_1_NRI$latitude
NRI_turnover_ratios$longitude <- combined_week_1_NRI$longitude
#add migratory proportion values
NRI_turnover_ratios$mig_prop_week_1 <- combined_week_1_NRI$proportion
NRI_turnover_ratios$mig_prop_week_26 <- combined_week_26_NRI$proportion
#calculate turnover in migratory proportion
NRI_turnover_ratios$mig_prop_turnover <- NRI_turnover_ratios$mig_prop_week_26 - NRI_turnover_ratios$mig_prop_week_1
#add individual community identifier
NRI_turnover_ratios$community <- combined_week_1_NRI$community

#export ecoregional seasonal turnover processed data
write.csv(NRI_turnover_ratios, "output/tables/NRI_turnover_ratios.csv")
write.csv(NTI_turnover_ratios, "output/tables/NTI_turnover_ratios.csv")

#read back in if already processed
NRI_turnover_ratios <- read.csv("output/tables/NRI_turnover_ratios.csv", row.names=1)
NTI_turnover_ratios <- read.csv("output/tables/NTI_turnover_ratios.csv", row.names=1)

##### CALCULATE SEASONAL TURNOVER VALUES IN PHYLOGENETIC COMMUNITY METRICS FOR CONTINENTAL COMMUNITIES #####

#GO BACK: get rid of those redo values
#reread in NRI and NTI ecoregional data generated from weeks 1 and weeks 26
NRI_week_26_continental <- read.csv("output/tables/26/passerine.ses.mpd.26.csv")
NTI_week_26_continental <- read.csv("output/tables/26/passerine.ses.mntd.26.csv")
NRI_week_1_continental <- read.csv("output/tables/1/passerine.ses.mpd.1_redo.csv")
NTI_week_1_continental <- read.csv("output/tables/1/passerine.ses.mntd.1_redo.csv")

#rename community column as "community" if it has been renamed as "X" (skip if not necessary)
NTI_week_1_continental <- NTI_week_1_continental %>%
  rename(community = X)
NRI_week_1_continental <- NRI_week_1_continental %>%
  rename(community = X)
NTI_week_26_continental <- NTI_week_26_continental %>%
  rename(community = X)
NRI_week_26_continental <- NRI_week_26_continental %>%
  rename(community = X)

### in case there are more values in week 1 than week 26, filter such that there are shared communities only
NTI_week_1_continental <- NTI_week_1_continental %>%
  filter(community %in% NTI_week_26_continental$community)
NRI_week_1_continental <- NRI_week_1_continental %>%
  filter(community %in% NRI_week_26_continental$community)

### PROCESS NTI DATA

#create data frames representing weeks 1 and 26 to initalize seasonal turnover dataframe
NTI_week_1_continental_df <- rep(NA, nrow(NTI_week_1_continental))
NTI_week_26_continental_df <- rep(NA, nrow(NTI_week_26_continental))
#merge into seasonal turnover into one dataframe
NTI_turnover_ratios_continental <- data.frame(NTI_week_1_continental_df, NTI_week_26_continental_df)
#add fields for species richness
NTI_turnover_ratios_continental$richness_week_1 <- rep(NA, nrow(NTI_week_1_continental))
NTI_turnover_ratios_continental$richness_week_26 <- rep(NA, nrow(NTI_week_1_continental))
#add field for seasonal turnover in NTI
NTI_turnover_ratios_continental$NTI_turnover <- rep(NA, nrow(NTI_week_1_continental))
#add field for seasonal turnover in species richness
NTI_turnover_ratios_continental$richness_turnover <- rep(NA, nrow(NTI_week_1_continental))
#add fields for BCR and geographic coordinates
NTI_turnover_ratios_continental$BCR <- rep(NA, nrow(NTI_week_1_continental))
NTI_turnover_ratios_continental$latitude <- rep(NA, nrow(NTI_week_1_continental))
NTI_turnover_ratios_continental$longitude <- rep(NA, nrow(NTI_week_1_continental))
#add fields for migratory proportion
NTI_turnover_ratios_continental$mig_prop_week_1 <- rep(NA, nrow(NTI_week_1_continental))
NTI_turnover_ratios_continental$mig_prop_week_26 <- rep(NA, nrow(NTI_week_1_continental))
#add field for seasonal turnover in migratory proportion
NTI_turnover_ratios_continental$mig_prop_turnover <- rep(NA, nrow(NTI_week_1_continental))
#add field for individual community identifiers
NTI_turnover_ratios_continental$community <- rep(NA, nrow(NTI_week_1_continental))

#populate dataframe generated above

#add NTI values
NTI_turnover_ratios_continental$NTI_week_1_continental <- NTI_week_1_continental$NTI
NTI_turnover_ratios_continental$NTI_week_26_continental <- NTI_week_26_continental$NTI
#add species richness values
NTI_turnover_ratios_continental$richness_week_1 <- NTI_week_1_continental$ntaxa
NTI_turnover_ratios_continental$richness_week_26 <- NTI_week_26_continental$ntaxa
#calculate turnover in NTI and species richness
NTI_turnover_ratios_continental$NTI_turnover <- NTI_turnover_ratios_continental$NTI_week_26_continental - NTI_turnover_ratios_continental$NTI_week_1_continental
NTI_turnover_ratios_continental$richness_turnover <- NTI_turnover_ratios_continental$richness_week_26 - NTI_turnover_ratios_continental$richness_week_1
#populate BCR and geographic coordinates
NTI_turnover_ratios_continental$BCR <- NTI_week_1_continental$BCR
NTI_turnover_ratios_continental$latitude <- NTI_week_1_continental$latitude
NTI_turnover_ratios_continental$longitude <- NTI_week_1_continental$longitude
#add migratory proportion values
NTI_turnover_ratios_continental$mig_prop_week_1 <- NTI_week_1_continental$proportion
NTI_turnover_ratios_continental$mig_prop_week_26 <- NTI_week_26_continental$proportion
#calculate turnover in migratory proportion
NTI_turnover_ratios_continental$mig_prop_turnover <- NTI_turnover_ratios_continental$mig_prop_week_26 - NTI_turnover_ratios_continental$mig_prop_week_1
#add individual community identifier
NTI_turnover_ratios_continental$community <- NTI_week_1_continental$community


### PROCESS NRI DATA

#create data frames representing weeks 1 and 26 to initalize seasonal turnover dataframe
NRI_week_1_continental_df <- rep(NA, nrow(NRI_week_1_continental))
NRI_week_26_continental_df <- rep(NA, nrow(NRI_week_26_continental))
#merge into seasonal turnover into one dataframe
NRI_turnover_ratios_continental <- data.frame(NRI_week_1_continental_df, NRI_week_26_continental_df)
#add fields for species richness
NRI_turnover_ratios_continental$richness_week_1 <- rep(NA, nrow(NRI_week_1_continental))
NRI_turnover_ratios_continental$richness_week_26 <- rep(NA, nrow(NRI_week_1_continental))
#add field for seasonal turnover in NRI
NRI_turnover_ratios_continental$NRI_turnover <- rep(NA, nrow(NRI_week_1_continental))
#add field for seasonal turnover in species richness
NRI_turnover_ratios_continental$richness_turnover <- rep(NA, nrow(NRI_week_1_continental))
#add fields for BCR and geographic coordinates
NRI_turnover_ratios_continental$BCR <- rep(NA, nrow(NRI_week_1_continental))
NRI_turnover_ratios_continental$latitude <- rep(NA, nrow(NRI_week_1_continental))
NRI_turnover_ratios_continental$longitude <- rep(NA, nrow(NRI_week_1_continental))
#add fields for migratory proportion
NRI_turnover_ratios_continental$mig_prop_week_1 <- rep(NA, nrow(NRI_week_1_continental))
NRI_turnover_ratios_continental$mig_prop_week_26 <- rep(NA, nrow(NRI_week_1_continental))
#add field for seasonal turnover in migratory proportion
NRI_turnover_ratios_continental$mig_prop_turnover <- rep(NA, nrow(NRI_week_1_continental))
#add field for individual community identifiers
NRI_turnover_ratios_continental$community <- rep(NA, nrow(NRI_week_1_continental))

#populate dataframe generated above

#add NRI values
NRI_turnover_ratios_continental$NRI_week_1_continental <- NRI_week_1_continental$NRI
NRI_turnover_ratios_continental$NRI_week_26_continental <- NRI_week_26_continental$NRI
#add species richness values
NRI_turnover_ratios_continental$richness_week_1 <- NRI_week_1_continental$ntaxa
NRI_turnover_ratios_continental$richness_week_26 <- NRI_week_26_continental_df$ntaxa
#calculate turnover in NRI and species richness
NRI_turnover_ratios_continental$NRI_turnover <- NRI_turnover_ratios_continental$NRI_week_26_continental - NRI_turnover_ratios_continental$NRI_week_1_continental
NRI_turnover_ratios_continental$richness_turnover <- NRI_turnover_ratios_continental$richness_week_26 - NRI_turnover_ratios_continental$richness_week_1
#populate BCR and geographic coordinates
NRI_turnover_ratios_continental$BCR <- NRI_week_1_continental$BCR
NRI_turnover_ratios_continental$latitude <- NRI_week_1_continental$latitude
NRI_turnover_ratios_continental$longitude <- NRI_week_1_continental$longitude
#add migratory proportion values
NRI_turnover_ratios_continental$mig_prop_week_1 <- NRI_week_1_continental$proportion
NRI_turnover_ratios_continental$mig_prop_week_26 <- NRI_week_26_continental_df$proportion
#calculate turnover in migratory proportion
NRI_turnover_ratios_continental$mig_prop_turnover <- NRI_turnover_ratios_continental$mig_prop_week_26 - NRI_turnover_ratios_continental$mig_prop_week_1
#add individual community identifier
NRI_turnover_ratios_continental$community <- NRI_week_1_continental$community

#export continental seasonal turnover processed data
write.csv(NRI_turnover_ratios_continental, "output/tables/NRI_turnover_ratios_continental.csv")
write.csv(NTI_turnover_ratios_continental, "output/tables/NTI_turnover_ratios_continental.csv")

#read back in if already processed
NRI_turnover_ratios_continental <- read.csv("output/tables/NRI_turnover_ratios_continental.csv", row.names=1)
NTI_turnover_ratios_continental <- read.csv("output/tables/NTI_turnover_ratios_continental.csv", row.names=1)

##### GENERATE PLOTS FOR SEASONAL TURNOVER IN COMMUNITY PHYLOGENETIC METRICS VS. SPECIES RICHNESS AND MIGRATORY PROPORTION (FIGURE 4) #####

#plot NTI turnover vs. richness turnover
NTI_richness_turnover_plot_continental <- ggplot(data = NTI_turnover_ratios_continental, aes(x=richness_turnover, y=NTI_turnover)) + 
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNTI vs. Richness Differences (June - January)") + 
  labs(x = "Richness Difference", y = "NTI Difference") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE)
ggMarginal(NTI_richness_turnover_plot_continental, fill = "green")

#plot NRI turnover vs. richness turnover
NRI_richness_turnover_plot_continental <- ggplot(data = NRI_turnover_ratios_continental, aes(x=richness_turnover, y=NRI_turnover)) +
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNRI vs. Richness Differences (June - January)") + 
  labs(x = "Richness Difference", y = "NRI Difference") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE)
ggMarginal(NRI_richness_turnover_plot_continental, fill = "green")

#plot NTI turnover vs. mig_prop turnover
NTI_mig_prop_turnover_plot_continental <- ggplot(data = NTI_turnover_ratios_continental, aes(x=mig_prop_turnover, y=NTI_turnover)) +
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNTI vs. Proportion Differences (June - January)") + 
  labs(x = "Proportion Difference", y = "NTI Difference") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE)
ggMarginal(NTI_mig_prop_turnover_plot_continental, fill = "green")

#plot NRI turnover vs. mig_prop turnover
NRI_mig_prop_turnover_plot_continental <- ggplot(data = NRI_turnover_ratios_continental, aes(x=mig_prop_turnover, y=NRI_turnover)) +
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNRI vs. Proportion Differences (June - January)") + 
  labs(x = "Proportion Difference", y = "NRI Difference") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE)
ggMarginal(NRI_mig_prop_turnover_plot_continental, fill = "green")


### BCR
#plot NTI turnover vs. richness turnover
NTI_richness_turnover_plot <- ggplot(data = NTI_turnover_ratios, aes(x=richness_turnover, y=NTI_turnover)) +
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNTI vs. Richness Differences (June - January)") + 
  labs(x = "Richness Difference", y = "NTI Difference") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE)
ggMarginal(NTI_richness_turnover_plot, fill = "green")

#plot NRI turnover vs. richness turnover
NRI_richness_turnover_plot <- ggplot(data = NRI_turnover_ratios, aes(x=richness_turnover, y=NRI_turnover)) +
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNRI vs. Richness Differences (June - January)") + 
  labs(x = "Richness Difference", y = "NRI Difference") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE)
ggMarginal(NRI_richness_turnover_plot, fill = "green")

#plot NTI turnover vs. mig_prop turnover
NTI_mig_prop_turnover_plot <- ggplot(data = NTI_turnover_ratios, aes(x=mig_prop_turnover, y=NTI_turnover)) +
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNTI vs. Proportion Differences (June - January)") + 
  labs(x = "Proportion Difference", y = "NTI Difference") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE)
ggMarginal(NTI_mig_prop_turnover_plot, fill = "green")

#plot NRI turnover vs. mig_prop turnover
NRI_mig_prop_turnover_plot <- ggplot(data = NRI_turnover_ratios, aes(x=mig_prop_turnover, y=NRI_turnover)) +
  geom_point(aes(color = latitude), size=4) + 
  scale_colour_viridis(name = "Latitude", direction = 1) + 
  ggtitle("\nNRI vs. Proportion Differences (June - January)") + 
  labs(x = "Proportion Difference", y = "NRI Difference") +
  theme(text = element_text(size = 20), legend.key.size = unit(1, 'cm')) +
  geom_smooth(method = "lm", se=TRUE)
ggMarginal(NRI_mig_prop_turnover_plot, fill = "green")


### try another new method: sdmTMB

#test out sdmTMB
head(NTI_turnover_ratios)

mesh <- make_mesh(NTI_turnover_ratios, xy_cols = c("longitude", "latitude"), cutoff = 50)
plot(mesh)

### regional
### regional NTI
#NTI ~ richness
sdmTMB_NTI_richness <- sdmTMB(
  formula = NTI_turnover ~ richness_turnover,
  data = NTI_turnover_ratios,
  mesh = mesh,
  spatial = "on"  # enables spatial random field
)
tidy(sdmTMB_NTI_richness, conf.int = TRUE)

#NTI ~ mig. prop.
sdmTMB_NTI_mig_prop <- sdmTMB(
  formula = NTI_turnover ~ mig_prop_turnover,
  data = NTI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NTI_mig_prop, conf.int = TRUE)

#NTI ~ richness + mig. prop.
sdmTMB_NTI_mig_prop_richness <- sdmTMB(
  formula = NTI_turnover ~ mig_prop_turnover + richness_turnover,
  data = NTI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NTI_mig_prop_richness, conf.int = TRUE)

#NTI ~ richness + mig. prop. + interaction
sdmTMB_NTI_mig_prop_richness_interaction <- sdmTMB(
  formula = NTI_turnover ~ mig_prop_turnover + richness_turnover + mig_prop_turnover*richness_turnover,
  data = NTI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NTI_mig_prop_richness_interaction, conf.int = TRUE)

#NTI ~ null
sdmTMB_NTI_null <- sdmTMB(
  formula = NTI_turnover ~ 1,
  data = NTI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NTI_null, conf.int = TRUE)

#create function to compare sdmTMB models
compare_sdmTMB_models <- function(..., sort = TRUE) {
  models <- list(...)
  model_names <- as.character(substitute(list(...)))[-1]
  
  out <- data.frame(
    Model = model_names,
    AIC = sapply(models, AIC),
    BIC = sapply(models, BIC),
    logLik = sapply(models, logLik)
  )
  
  out$Delta_AIC <- out$AIC - min(out$AIC)
  out$Weight <- exp(-0.5 * out$Delta_AIC)
  out$Weight <- out$Weight / sum(out$Weight)
  
  if (sort) out <- out[order(out$AIC), ]
  rownames(out) <- NULL
  return(out)
}

#model selection
compare_sdmTMB_models(sdmTMB_NTI_richness, sdmTMB_NTI_mig_prop, sdmTMB_NTI_mig_prop_richness,
                      sdmTMB_NTI_mig_prop_richness_interaction, sdmTMB_NTI_null)


###regional NRI
### regional NRI
#NRI ~ richness
sdmTMB_NRI_richness <- sdmTMB(
  formula = NRI_turnover ~ richness_turnover,
  data = NRI_turnover_ratios,
  mesh = mesh,
  spatial = "on"  # enables spatial random field
)
tidy(sdmTMB_NRI_richness, conf.int = TRUE)

#NRI ~ mig. prop.
sdmTMB_NRI_mig_prop <- sdmTMB(
  formula = NRI_turnover ~ mig_prop_turnover,
  data = NRI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NRI_mig_prop, conf.int = TRUE)

#NRI ~ richness + mig. prop.
sdmTMB_NRI_mig_prop_richness <- sdmTMB(
  formula = NRI_turnover ~ mig_prop_turnover + richness_turnover,
  data = NRI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NRI_mig_prop_richness, conf.int = TRUE)

#NRI ~ richness + mig. prop. + interaction
sdmTMB_NRI_mig_prop_richness_interaction <- sdmTMB(
  formula = NRI_turnover ~ mig_prop_turnover + richness_turnover + mig_prop_turnover*richness_turnover,
  data = NRI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NRI_mig_prop_richness_interaction, conf.int = TRUE)

#NRI ~ null
sdmTMB_NRI_null <- sdmTMB(
  formula = NRI_turnover ~ 1,
  data = NRI_turnover_ratios,
  mesh = mesh,
  spatial = "on"
)
tidy(sdmTMB_NRI_null, conf.int = TRUE)

#model selection
compare_sdmTMB_models(sdmTMB_NRI_richness, sdmTMB_NRI_mig_prop, sdmTMB_NRI_mig_prop_richness,
                      sdmTMB_NRI_mig_prop_richness_interaction, sdmTMB_NRI_null)





### continental
### continental NTI
mesh_continental <- make_mesh(NTI_turnover_ratios_continental, xy_cols = c("latitude", "longitude"), cutoff = 50)

#NTI ~ richness
sdmTMB_NTI_richness_continental <- sdmTMB(
  formula = NTI_turnover ~ richness_turnover,
  data = NTI_turnover_ratios_continental,
  mesh = mesh_continental,
  spatial = "on"  # enables spatial random field
)
tidy(sdmTMB_NTI_richness_continental, conf.int = TRUE)

#NTI ~ mig. prop.
sdmTMB_NTI_mig_prop_continental <- sdmTMB(
  formula = NTI_turnover ~ mig_prop_turnover,
  data = NTI_turnover_ratios_continental,
  mesh = mesh_continental,
  spatial = "on"
)
tidy(sdmTMB_NTI_mig_prop_continental, conf.int = TRUE)

#NTI ~ richness + mig. prop.
sdmTMB_NTI_mig_prop_richness_continental <- sdmTMB(
  formula = NTI_turnover ~ mig_prop_turnover + richness_turnover,
  data = NTI_turnover_ratios_continental,
  mesh = mesh_continental,
  spatial = "on"
)
tidy(sdmTMB_NTI_mig_prop_richness_continental, conf.int = TRUE)

#NTI ~ richness + mig. prop. + interaction
sdmTMB_NTI_mig_prop_richness_interaction_continental <- sdmTMB(
  formula = NTI_turnover ~ mig_prop_turnover + richness_turnover + mig_prop_turnover*richness_turnover,
  data = NTI_turnover_ratios_continental,
  mesh = mesh_continental,
  spatial = "on"
)
tidy(sdmTMB_NTI_mig_prop_richness_interaction_continental, conf.int = TRUE)

#NTI ~ null
sdmTMB_NTI_null_continental <- sdmTMB(
  formula = NTI_turnover ~ 1,
  data = NTI_turnover_ratios_continental,
  mesh = mesh_continental,
  spatial = "on"
)
tidy(sdmTMB_NTI_null_continental, conf.int = TRUE)

#model selection
compare_sdmTMB_models(sdmTMB_NTI_richness_continental, sdmTMB_NTI_mig_prop_continental, sdmTMB_NTI_mig_prop_richness_continental,
                      sdmTMB_NTI_mig_prop_richness_interaction_continental, sdmTMB_NTI_null_continental)



### continental NRI

#NRI ~ richness
sdmTMB_NRI_richness_continental <- sdmTMB(
  formula = NRI_turnover ~ richness_turnover,
  data = NRI_turnover_ratios_continental,
  mesh = mesh_continental,
  spatial = "on"  # enables spatial random field
)
tidy(sdmTMB_NRI_richness_continental, conf.int = TRUE)

#NRI ~ mig. prop.
sdmTMB_NRI_mig_prop_continental <- sdmTMB(
  formula = NRI_turnover ~ mig_prop_turnover,
  data = NRI_turnover_ratios_continental,
  mesh = mesh_continental,
  spatial = "on"
)
tidy(sdmTMB_NRI_mig_prop_continental, conf.int = TRUE)

#NRI ~ richness + mig. prop.
sdmTMB_NRI_mig_prop_richness_continental <- sdmTMB(
  formula = NRI_turnover ~ mig_prop_turnover + richness_turnover,
  data = NRI_turnover_ratios_continental,
  mesh = mesh_continental,
  spatial = "on"
)
tidy(sdmTMB_NRI_mig_prop_richness_continental, conf.int = TRUE)

#NRI ~ richness + mig. prop. + interaction
sdmTMB_NRI_mig_prop_richness_interaction_continental <- sdmTMB(
  formula = NRI_turnover ~ mig_prop_turnover + richness_turnover + mig_prop_turnover*richness_turnover,
  data = NRI_turnover_ratios_continental,
  mesh = mesh_continental,
  spatial = "on"
)
tidy(sdmTMB_NRI_mig_prop_richness_interaction_continental, conf.int = TRUE)

#NRI ~ null
sdmTMB_NRI_null_continental <- sdmTMB(
  formula = NRI_turnover ~ 1,
  data = NRI_turnover_ratios_continental,
  mesh = mesh_continental,
  spatial = "on"
)
tidy(sdmTMB_NRI_null_continental, conf.int = TRUE)

#model selection
compare_sdmTMB_models(sdmTMB_NRI_richness_continental, sdmTMB_NRI_mig_prop_continental, sdmTMB_NRI_mig_prop_richness_continental,
                      sdmTMB_NRI_mig_prop_richness_interaction_continental, sdmTMB_NRI_null_continental)

#test for spatial autocorrelation

#NTI turnover ratios
test_NTI_turnover_ratios <- NTI_turnover_ratios[,-c(3,9)]
test_NTI_turnover_ratios <- test_NTI_turnover_ratios[complete.cases(test_NTI_turnover_ratios),]
coords <- cbind(test_NTI_turnover_ratios$longitude, test_NTI_turnover_ratios$latitude)
nb <- dnearneigh(coords, 0, 5)
lw <- nb2listw(nb)
model <- gls(NTI_turnover ~ richness_turnover, test_NTI_turnover_ratios, method = "REML")
moran.test(residuals(model), lw)
#Moran's I statistic standard deviate: 509.48
#p-value: <2.2 e -16
#Moran I statistic: 1.3426 e -01

#NRI turnover ratios
test_NRI_turnover_ratios <- NRI_turnover_ratios[,-c(3,9)]
test_NRI_turnover_ratios <- test_NRI_turnover_ratios[complete.cases(test_NRI_turnover_ratios),]
coords2 <- cbind(test_NRI_turnover_ratios$longitude, test_NRI_turnover_ratios$latitude)
nb2 <- dnearneigh(coords2, 0, 5)
lw2<- nb2listw(nb2)
model2 <- gls(NRI_turnover ~ richness_turnover, test_NRI_turnover_ratios, method = "REML")
moran.test(residuals(model2), lw2)
#Moran's I statistic standard deviate: 933.28
#p-value: <2.2 e -16
#Moran I statistic: 2.4597 e -01

#NTI turnover ratios continental
test_NTI_turnover_ratios_continental <- NTI_turnover_ratios_continental[,-c(3,9)]
test_NTI_turnover_ratios_continental <- test_NTI_turnover_ratios_continental[complete.cases(test_NTI_turnover_ratios_continental),]
coords3 <- cbind(test_NTI_turnover_ratios_continental$longitude, test_NTI_turnover_ratios_continental$latitude)
nb3 <- dnearneigh(coords3, 0, 5)
lw3<- nb2listw(nb3)
model3 <- gls(NTI_turnover ~ richness_turnover, test_NTI_turnover_ratios_continental, method = "REML")
moran.test(residuals(model3), lw3)
#Moran's I statistic standard deviate: 753.83
#p-value: <2.2 e -16
#Moran I statistic: 1.98667 e -01

#NTI turnover ratios continental
test_NRI_turnover_ratios_continental <- NRI_turnover_ratios_continental[,-c(3,9)]
test_NRI_turnover_ratios_continental <- test_NRI_turnover_ratios_continental[complete.cases(test_NRI_turnover_ratios_continental),]
coords4 <- cbind(test_NRI_turnover_ratios_continental$longitude, test_NRI_turnover_ratios_continental$latitude)
nb4 <- dnearneigh(coords4, 0, 5)
lw4<- nb2listw(nb4)
model4 <- gls(NRI_turnover ~ richness_turnover, test_NRI_turnover_ratios_continental, method = "REML")
moran.test(residuals(model4), lw4)
#Moran's I statistic standard deviate: 1199.1
#p-value: <2.2 e -16
#Moran I statistic: 3.1606 e -01


### Consider re-doing Figure 4

crs_102004 <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

na_basemap <- ne_countries(continent = "North America", returnclass = "sf")
na_basemap_proj <- st_transform(na_basemap, crs = crs_102004)

NRI_turnover_ratios_sf <- st_as_sf(NRI_turnover_ratios, coords = c("longitude", "latitude"), crs = 4326)
NRI_turnover_ratios_proj <- st_transform(NRI_turnover_ratios_sf, crs = crs_102004)
NRI_turnover_ratios_proj <- NRI_turnover_ratios_proj %>%
  mutate(X = st_coordinates(.)[,1],
         Y = st_coordinates(.)[,2])

#plot migratory proportion difference
ggplot() +
  geom_sf(data = na_basemap_proj, fill = "gray95", color = "gray60") +
  geom_sf(data = NRI_turnover_ratios_proj, aes(color = mig_prop_turnover), size = 0.75) +
  coord_sf(crs = crs_102004,
           xlim = c(-3550000, 3000000),
           ylim = c(-2200000, 3000000)) +
  ggtitle("Seasonal Migratory Proportion Difference (June - January)") +
  theme(text = element_text(size = 20), legend.title = element_text(hjust = 0.5)) +
  scale_color_gradient2("Mig. Prop. \nDiff.",
                        midpoint = 0,
                        low = "red",
                        mid = "white",
                        high = "blue")

#plot species richness difference
ggplot() +
  geom_sf(data = na_basemap_proj, fill = "gray95", color = "gray60") +
  geom_sf(data = NRI_turnover_ratios_proj, aes(color = richness_turnover), size = 0.75) +
  coord_sf(crs = crs_102004,
           xlim = c(-3550000, 3000000),
           ylim = c(-2200000, 3000000)) +
  ggtitle("Seasonal Species Richness Difference (June - January)") +
  theme(text = element_text(size = 20), legend.title = element_text(hjust = 0.5)) +
  scale_color_gradient2("Richness \nDiff.",
                        midpoint = 0,
                        low = "red",
                        mid = "white",
                        high = "blue") 

#plot regional NRI turnover
ggplot() +
  geom_sf(data = na_basemap_proj, fill = "gray95", color = "gray60") +
  geom_sf(data = NRI_turnover_ratios_proj, aes(color = NRI_turnover), size = 0.75) +
  coord_sf(crs = crs_102004,
           xlim = c(-3550000, 3000000),
           ylim = c(-2200000, 3000000)) +
  ggtitle("Seasonal NRI Difference (June - January), Regional Pool") +
  theme(text = element_text(size = 20), legend.title = element_text(hjust = 0.5)) +
  scale_color_gradient2("NRI Diff.",
                        midpoint = 0,
                        low = "red",
                        mid = "white",
                        high = "blue")

#add NTI ratios to this projection dataframe
NRI_turnover_ratios_proj$NTI_turnover <- NTI_turnover_ratios$NTI_turnover

#plot regional NTI turnover
ggplot() +
  geom_sf(data = na_basemap_proj, fill = "gray95", color = "gray60") +
  geom_sf(data = NRI_turnover_ratios_proj, aes(color = NTI_turnover), size = 0.75) +
  coord_sf(crs = crs_102004,
           xlim = c(-3550000, 3000000),
           ylim = c(-2200000, 3000000)) +
  ggtitle("Seasonal NTI Difference (June - January), Regional Pool") +
  theme(text = element_text(size = 20), legend.title = element_text(hjust = 0.5)) +
  scale_color_gradient2("NTI Diff.",
                        midpoint = 0,
                        low = "red",
                        mid = "white",
                        high = "blue")

#continental data plotting
NRI_turnover_ratios_continental_sf <- st_as_sf(NRI_turnover_ratios_continental, coords = c("longitude", "latitude"), crs = 4326)
NRI_turnover_ratios_continental_proj <- st_transform(NRI_turnover_ratios_continental_sf, crs = crs_102004)
NRI_turnover_ratios_continental_proj <- NRI_turnover_ratios_continental_proj %>%
  mutate(X = st_coordinates(.)[,1],
         Y = st_coordinates(.)[,2])


#plot continental NRI turnover
ggplot() +
  geom_sf(data = na_basemap_proj, fill = "gray95", color = "gray60") +
  geom_sf(data = NRI_turnover_ratios_continental_proj, aes(color = NRI_turnover), size = 0.75) +
  coord_sf(crs = crs_102004,
           xlim = c(-3550000, 3000000),
           ylim = c(-2200000, 3000000)) +
  ggtitle("Seasonal NRI Difference (June - January), Continental Pool") +
  theme(text = element_text(size = 20), legend.title = element_text(hjust = 0.5)) +
  scale_color_gradient2("NRI Diff.",
                        midpoint = 0,
                        low = "red",
                        mid = "white",
                        high = "blue")

#plot continental NTI turnover

NRI_turnover_ratios_continental_proj$NTI_turnover <- NTI_turnover_ratios_continental$NTI_turnover

ggplot() +
  geom_sf(data = na_basemap_proj, fill = "gray95", color = "gray60") +
  geom_sf(data = NRI_turnover_ratios_continental_proj, aes(color = NTI_turnover), size = 0.75) +
  coord_sf(crs = crs_102004,
           xlim = c(-3550000, 3000000),
           ylim = c(-2200000, 3000000)) +
  ggtitle("Seasonal NTI Difference (June - January), Continental Pool") +
  theme(text = element_text(size = 20), legend.title = element_text(hjust = 0.5)) +
  scale_color_gradient2("NTI Diff.",
                        midpoint = 0,
                        low = "red",
                        mid = "white",
                        high = "blue")
