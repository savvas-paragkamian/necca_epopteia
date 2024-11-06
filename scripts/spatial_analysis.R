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

## greece 
##
## borders for maps
#world <- ne_countries(scale = "medium", returnclass = "sf")

#country_name <- "Greece"
#greece <- world[world$name == country_name, ]

hellenic_borders_shp <- sf::st_read("../spatial_data/hellenic_borders/hellenic_borders.shp")

#### World Clim , Bioclim Variables
###
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
world_clim_directory <- "/Users/talos/Documents/spatial_data/world_clim/wc2.1_30s_bio/"
output_directory <- "/Users/talos/Documents/programming_projects/necca_epopteia/spatial_data/world_clim_greece/"

world_clim_files <- list.files(world_clim_directory)

bbox_polygon <- st_as_sf(st_as_sfc(st_bbox(hellenic_borders_shp)))

for (f in world_clim_files) {
    
    if (grepl("*.tif$", f)) {
        
        #read_raster
        path_raster <- paste0(world_clim_directory,f,sep="")
        raster_tmp <- rast(path_raster)
        
        crete_raster <- terra::crop(raster_tmp, bbox_polygon)
        crete_raster <- terra::project(crete_raster, wgs84)
        output_raster <- paste0(output_directory, "crete_",f,sep="")
        print(output_raster)
        terra::writeRaster(crete_raster, output_raster,overwrite=TRUE)

        rm(path_raster,raster_tmp,crete_raster,output_raster)

    }else{
        
        print(f, " not a tif")
        next
    }
}
