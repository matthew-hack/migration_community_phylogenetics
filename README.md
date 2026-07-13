# Migration community phylogenetics
Pipeline for data generation and analysis of community phylogenetics of migration, using eBird Status and Trends data.

There are six scripts associated with this text:

-- 1_SDM_download_and_community_formation.R <br />
This script accesses the eBird Status and Trends data package in order to download the eBird Status and Trends data for each species included in the dataset. The script proceeds to partition the eBirdst data into weekly distributional rasters, which are then aggregated into community matrices for each week of the year depicting the abundance of each individual species.

-- 2_weekly_community_passerine_phylogenetic_metrics.R <br />
This script runs all of the community phylogenetic metric analyses for one week of the year (the first week) -- the same script can be utilized for the other 51 weeks of the year. NRI and NTI metrics are generated for each community at two species pools (a continental and ecoregional pool). The script also creates visualizatoins for the geospatial distributions of NRI and NTI metrics at both species pools.

-- 3_community_migratory_proportion_calculations.R <br />
This script defines a function that counts the number of migratory species in each community during each week of the year, and then calculates the proportion of migratory species in each community. The script gives an example of the pipeline for a single week (week 1), which can be substituted for the remaining weeks of the year. The script proceeds to apply the migratory proportion data to the previously-generated community phylogenetic metric output in the previous script.

-- 4_exploratory_values_and_figures.R <br />
This script generates the values used in the Results section of the manuscript, along with several exploratory figures depicting the range and geographic distribution of NRI and NTI values, along with variation in species richness and migratory proportion.

-- 5_manuscript_data_analysis.R <br />
This script processes seasonal turnover (from winter to summer) in community phylogenetic metrics, along with species richness and migratory proportion, across both species pool levels. Spatial GLMMs are carried out to analyze the seasonality of community phylogenetic metrics through space, and to model seasonal turnover in community phylogenetic metrics against other factors. The script also analyzes the phylogenetic distance of NRI and NTI values within summer and winter (SF1) and generates figures depicting seasonal turnover in phylogenetic community metrics.

-- 6_seasonality_analysis.R <br />
This script downloads and analyzes bioclimatic variables, and tests the influence of annual temperature range and seasonality in shaping seasonal variation in community phylogenetic metrics.
