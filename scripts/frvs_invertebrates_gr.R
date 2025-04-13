#!/usr/bin/Rscript

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


########################## Flowchart for FRVs ##########################
########################## Parnassius apollo ###########################


parnassius_dist <- sf::st_read("../data/Parnassius apollo AP 2019/AP_Papollo_Distribution_LAEA.shp") |>
    st_transform(crs="WGS84")

species_samples_art17_parnasious <- species_samples_art17_sf |>
    filter(species=="Parnassius apollo") 

gbif_parnasious_gr <- st_intersection(gbif_parnasious,greece_regions)


