#!/usr/bin/Rscript

## Script name: db_standardization.R
##
## Purpose of script:
## Find the differences in the structure and content of 
## the different datasets of invertebrate data.
## More specifically: NECCA db, ENVECO db, Monitoring data
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-12-15

library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(units)

###
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# NECCA gtb
# Define the directory containing GTB files
dir_path <- "../data/necca_internal_extract.gdb"

# layers of the database
gtb_layers <- st_layers(dir_path)

# Read all GTB files into a list of sf objects
gtb_data_list <- purrr::map(gtb_layers$name,~st_read(dsn=dir_path,layer=.))

# Optionally, combine them into a single sf object if structures are the same
## do not run. PROTOCOL_1A, GENIKES_INFO, PROTOCOL_1A__ATTACH, GENIKES_INFO__ATTACH
## don't have CRS
#gtb_combined <- do.call(rbind, gtb_data_list)
