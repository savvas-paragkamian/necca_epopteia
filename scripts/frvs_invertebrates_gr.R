#!/usr/bin/env Rscript

## Script name: dfrvs_invertebrates_gr.R
##
## Purpose of script: Calculate the Favourite Reference Values
## of the invertebrates in Greece
##
## Author: Savvas Paragkamian
##
## Date Created: 2025-03-13

library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(dismo)
library(units)

############################# Load Spatial Data ########################
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")
#hellenic_borders_shp <- sf::st_read("../spatial_data/hellenic_borders/hellenic_borders.shp")

gr_1km <- sf::st_read("../spatial_data/eea_1km/gr_1km.shp") |>
    st_transform(., crs="WGS84")

########################### Load Species Data ###########################
### Species occurrences enriched ######
species_occurrences_invertebrates <- read_delim("../results/species_occurrences_invertebrates.tsv",delim="\t")
species_occurrences_spatial <- read_delim("../results/species_occurrences_spatial.tsv",delim="\t")
species_samples_art17 <- read_delim("../results/species_samples_art17.tsv", delim="\t")

classification_species_gbif <- read_delim("../results/classification_species_gbif.tsv", delim="\t")
iucn_art17_invert_all <- read_delim("../results/iucn_art17_invert_all.tsv", delim="\t")

sspecies_dist_national_rep <- sf::st_read("../spatial_data/National report_2013_2018_shp/GR_Art17_species_distribution.shp")
species_dist_national_rep_sens <- sf::st_read("../spatial_data/National report_2013_2018_shp/GR_Art17_species_distribution_sensitive.shp")


########################## Flowchart for FRVs ##########################
# what is known for each species
species_info <- species_samples_art17 |>
    distinct(species) |> 
    left_join(classification_species_gbif, by=c("species"="species")) |>
    left_join(iucn_art17_invert_all, by=c("species"="scientificName"))



########################## Parnassius apollo ###########################


parnassius_dist <- sf::st_read("../data/Parnassius apollo AP 2019/AP_Papollo_Distribution_LAEA.shp") |>
    st_transform(crs="WGS84")

species_samples_art17_parnasious <- species_samples_art17_sf |>
    filter(species=="Parnassius apollo") 



