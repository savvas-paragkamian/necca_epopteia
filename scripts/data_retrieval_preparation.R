#!/usr/bin/Rscript

## Script name: data_retrieval_preparation.R
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-11-06

library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(taxize)
library(rgbif)
library(units)
library(vegan)
library(rnaturalearth)


###
###
greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")

######################## GBIF ########################
# GBIF retrieve data for all arthropod species that have been assessed in IUCN
### NOT run takes time. 
#species_names <- as.character(invertebrates_all_species$SPECIES_NAME)
#gbif_species <- get_gbifid(species_names,ask=F)


#species_gbif_df <- data.frame(sci_name=species_names, gbifid=gbif_species)
#write_delim(species_gbif_df, "../results/species_gbif_invertebrates.tsv", delim="\t")
#### takes even more time!!!
#classification_species <- classification(species_gbif_df$gbifid.ids, db = 'gbif')

#classification_species_d <- do.call(rbind, classification_species) |>
#    rownames_to_column(var="gbif") |> 
#    mutate(gbif = gsub("\\.(.*)","", gbif)) |>
#    dplyr::select(-id) |>
#    distinct() |>
#    na.omit(gbif) |>
#    pivot_wider(id_cols=gbif, names_from=rank, values_from=name ) |>
#    mutate(gbif=as.numeric(gbif)) 
#
#write_delim(classification_species_d, "../results/classification_species_gbif.tsv", delim="\t")
#
# Resolve names
## gnr_datasources() |> filter(title=="GBIF Backbone Taxonomy") id=11
#gnr_species <- gnr_resolve(species_gbif_df$sci_name)
#gnr_species_gbif <- gnr_resolve(species_gbif_df$sci_name, data_source_ids=11)
#
#write_delim(gnr_species, "../results/gnr_species_names.tsv", delim="\t")
#
########################## species occurrences ##############################
# GBIF occurrences
## need to set GBIF credential in the .Renviron file

#gbif_taxon_keys <- as.numeric(na.omit(species_gbif_df$gbifid.ids))
#
## run once to request the download from the server
#occ_download(
#pred_in("taxonKey", gbif_taxon_keys), # important to use pred_in
#pred("hasCoordinate", TRUE),
#pred("hasGeospatialIssue", FALSE),
#format = "SIMPLE_CSV"
#)

# to check the status of the download
# This is for the species of annex II 
# occ_download_wait('0026745-241024112534372') 
# the key for the gbif download of 268 invertegrate species 
# is 0018673-241107131044228 

#gbif_species_occ <- occ_download_get('0018673-241107131044228') |>
#    occ_download_import()

#write_delim(gbif_species_occ, "../results/gbif_invertebrate_species_occ.tsv", delim="\t")

gbif_species_occ <- read_delim("../results/gbif_invertebrate_species_occ.tsv", delim="\t")

gbif_species_occ_sf <- gbif_species_occ |> 
    st_as_sf(coords=c("decimalLatitude","decimalLongitude"),
             remove=F,
             crs="WGS84") |>
    st_transform(4326)

greece_regions_bbox <- st_as_sf(st_as_sfc(st_bbox(greece_regions)))

gbif_species_world_map <- ggplot() +
    geom_sf(world, mapping=aes()) +
    geom_point(gbif_species_occ_sf,
            mapping=aes(x=decimalLongitude,
                        y=decimalLatitude,
                        color=species),
            size=1.8,
            alpha=0.8,
            show.legend=T) +
    geom_sf(greece_regions, mapping=aes()) +
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

#ggsave("../figures/gbif_species_world_map.png", 
#       plot=gbif_species_world_map, 
#       height = 40, 
#       width = 60,
#       dpi = 300, 
#       units="cm",
#       device="png")

## Greece only

# Define the bounding box coordinates
xmin <- 19.37359
ymin <- 34.80202
xmax <- 29.64306
ymax <- 41.7485

gbif_species_occ_gr <- gbif_species_occ_sf |>
    filter(decimalLongitude > 19.37359 & decimalLongitude<29.64306,
           decimalLatitude>34.80202 & decimalLatitude < 41.7485 )

gbif_species_gr_map <- ggplot() +
    geom_point(gbif_species_occ_gr,
            mapping=aes(x=decimalLongitude,
                        y=decimalLatitude,
                        color=species),
            size=1.8,
            alpha=0.8,
            show.legend=T) +
    geom_sf(greece_regions, mapping=aes()) +
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../figures/gbif_species_gr_map.png", 
       plot=gbif_species_gr_map, 
       height = 40, 
       width = 40,
       dpi = 300, 
       units="cm",
       device="png")

################################## world Clim ####################################

## function to crop tif rasters based on bbox of other shapefile
## and make it wgs84 <- "EPSG:4326" 

#raster_path <- eurodem_file
#bbox_sf <- bbox_polygon
#output_path <- eurodem_gr_d
crop_raster_by_shp <- function(raster_path, shp, output_path) {
    # Read the raster
    raster_tmp <- rast(raster_path)
    bbox_polygon <- st_as_sf(st_as_sfc(st_bbox(shp)))
    # Check CRS match
    raster_crs <- crs(raster_tmp, proj=TRUE)
    bbox_crs <- st_crs(bbox_sf)$wkt

    if (raster_crs != bbox_crs) {
        message("Reprojecting bbox to match raster CRS...")
        bbox_sf <- st_transform(bbox_sf, crs = st_crs(raster_crs))
    }
    # Convert bbox to SpatVector for terra cropping
    #bbox_vect <- vect(bbox_sf)
    #Crop raster
    cropped_raster <- crop(raster_tmp, bbox_sf)
    # project to WGS84
    wgs84 <- "EPSG:4326"
    if (crs(cropped_raster) != wgs84) {
        cropped_raster <- project(cropped_raster, wgs84)
    }
    mask_raster <- mask(cropped_raster, shp)

    output_path <- paste0(output_path,"crop_",basename(raster_path),sep="")
  
    # Save the cropped raster
    terra::writeRaster(mask_raster, output_path, overwrite = TRUE)
  
    # Clean up
    rm(raster_tmp, cropped_raster, mask_raster)
    gc()
  
    message("Saved cropped raster to: ", output_path)
}

####### World clim #########
world_clim_directory <- "/Users/talos/Documents/spatial_data/world_clim/wc2.1_30s_bio/"
output_directory <- "/Users/talos/Documents/programming_projects/necca_epopteia/spatial_data/world_clim_greece/"

world_clim_files <- list.files(world_clim_directory, pattern = "\\.tif$", full.names = TRUE)

for (f in world_clim_files) {
    
    if (grepl("*.tif$", f)) {
        
        #read_raster
        path_raster <- f
        crop_raster_by_shp(path_raster,greece_regions,output_directory)
        
        rm(path_raster,raster_tmp,crete_raster,output_raster)

    }else{
        
        print(f, " not a tif")
        next
    }
}

################ HILDA plus

hilda_directory <- "/Users/talos/Documents/spatial_data/hildap_vGLOB-1.0_geotiff_wgs84/hildap_GLOB-v1.0_lulc-states/"
output_directory <- "/Users/talos/Documents/programming_projects/necca_epopteia/spatial_data/hilda_greece/"

hilda_directory_files <- list.files(hilda_directory, pattern = "\\.tif$", full.names = TRUE)


for (f in hilda_directory_files) {
    
    if (grepl("*.tif$", f)) {
        
        #read_raster
        path_raster <- f
        crop_raster_by_shp(path_raster,greece_regions,output_directory)
        
        rm(path_raster,raster_tmp,crete_raster,output_raster)

    }else{
        
        print(f, " not a tif")
        next
    }
}
#### EU DEM
###


eurodem_file <- "/Users/talos/Documents/spatial_data/euro-dem-tif/data/eurodem/eurodem.tif"
eurodem_gr_d <- "../spatial_data/eurodem_gr"
crop_raster_by_bbox(eurodem_file,bbox_polygon,eurodem_gr_d)

eu_dem <- rast(eurodem_file)

crete_raster <- terra::crop(raster_tmp, bbox_polygon)
crete_raster <- terra::project(crete_raster, wgs84)
output_raster <- paste0(output_directory, "crete_",f,sep="")
print(output_raster)
terra::writeRaster(crete_raster, output_raster,overwrite=TRUE)




