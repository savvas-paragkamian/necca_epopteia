#!/usr/bin/Rscript

## Script name: spatial_analysis.R
##
## Purpose of script:
##
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-11-06

library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(units)


########################## Load Data ##########################
###
### Species occurrences enriched ######
species_occurrences_art17_invertebrates <- read_delim("../results/species_occurrences_art17_invertebrates.tsv",delim="\t")
### greece 

greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")

#EEA 1km
eea_1km <- sf::st_read("../spatial_data/eea_reference_grid/gr_1km.shp")
eea_10km <- sf::st_read("../spatial_data/eea_reference_grid/gr_10km.shp")

### World Clim , Bioclim Variables
###
world_clim_directory <- "../spatial_data/world_clim_greece/"
world_clim_files <- list.files(world_clim_directory, pattern = "\\.tif$", full.names = TRUE)

world_clim_list <- lapply(world_clim_files, rast)

### Natura2000

N2000_v32 <- sf::st_read("../spatial_data/N2000_spatial_GR_2021_12_09_v32/N2000_spatial_GR_2021_12_09_v32.shp")

#01_Χαρτογράφηση χερσαίων Τ.Ο._EKXA








######################
# Define bounding box as sf object for visual checks if needed
bbox <- st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = 4326)
bbox_sf <- st_as_sfc(bbox)

# Generate 4 random points inside the bounding box
set.seed(123)
points_inside <- st_as_sf(data.frame(
  id = 1:4,
  x = runif(4, xmin, xmax),
  y = runif(4, ymin, ymax)
), coords = c("x", "y"), crs = 4326)

# Generate 4 random points outside the bounding box
points_outside <- st_as_sf(data.frame(
  id = 5:8,
  x = c(runif(2, xmin - 10, xmin - 1), runif(2, xmax + 1, xmax + 10)),
  y = c(runif(2, ymin - 10, ymin - 1), runif(2, ymax + 1, ymax + 10))
), coords = c("x", "y"), crs = 4326)

# Combine the points into one sf object for easy handling
all_points <- rbind(points_inside, points_outside)



