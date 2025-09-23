#!/usr/bin/Rscript

## Script name: data_retrieval_preparation.R
##
## Author: Savvas Paragkamian
## Purpose of script:
## preparation of different types of data, spatial, species occurrences,
## taxonomy. The script is interactive, i.e not to executed all at once.
##
## Date Created: 2024-11-06

library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(taxize)
library(rgbif)
library(units)
#library(rnaturalearth)
source("necca_spatial_functions.R")


###
###
greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")
greece_regions_vec <- vect(greece_regions)

greece_regions_bbox <- st_bbox(greece_regions) |> st_as_sfc()
greece_regions_bbox_sf <- st_sf(geometry = greece_regions_bbox, crs = st_crs(greece_regions))


greece_regions_bbox_LAEA <- greece_regions_bbox_sf |> 
    st_transform(3035)
########################## species information ##############################
# the species of the previous epopteia 2015

invertebrates_92_43 <- readxl::read_excel("../data/invertebrates_92_43.xlsx", skip=2, col_names=F)

colnames(invertebrates_92_43) <- c("no","area","SPECIES_NAME","SPECIES_ID","ORDER","CLASS","PRIORITY","ANNEX_II","ANNEX_IV","ANNEX_V","KD","POPULATION_TREND","POPULATION_SIZE_UNIT","OCCURRENCE","SD")

species_names <- unique(invertebrates_92_43$SPECIES_NAME) 

#### this list contains species names from multiple versions of Monitoring, EPOPTEIA I, EPOPTEIA II. 

species_names_combined <- c(
  "Apatura metis",
  "Astacus astacus",
  "Austropotamobious torrentium",
  "Austropotamobius torrentium",
  "Bolbelasmus unicornis",
  "Buprestis splendens",
  "Callimorpha (Euplagia, Panaxia) quadripunctaria",
  "Callimorpha quadripunctaria",
  "Catopta thrips",
  "Paracossulus thrips",
  "Cerambyx cerdo",
  "Coenagrion ornatum",
  "Cordulegaster heros",
  "Dioszeghyana schmidtii",
  "Eriogaster catax",
  "Euphydryas (Eurodryas, Hypodryas) aurinia",
  "Euphydryas aurinia",
  "Euplagia quadripunctaria",
  "Hirudo medicinalis",
  "Hirudo verbana",
  "Hyles hippophaes",
  "Lindenia tetraphylla",
  "Lucanus cervus",
  "Lycaena dispar",
  "Maculinea arion",
  "Morimus asper funereus",
  "Morimus asper",
  "Morimus funereus",
  "Ophiogomphus cecilia",
  "Osmoderma lassallei",
  "Osmoderma barnabita",
  "Osmoderma eremita",
  "Osmoderma eremita Complex",
  "Papilio alexanor",
  "Paracaloptenus caloptenoides",
  "Parnassius apollo",
  "Parnassius mnemosyne",
  "Polyommatus eroides",
  "Probaticus subrugosus",
  "Proserpinus proserpina",
  "Pseudophilotes bavius",
  "Rhysodes sulcatus",
  "Rosalia alpina",
  "Stenobothrus (Stenobothrodes) eurasius",
  "Stenobothrus eurasius",
  "Stylurus flavipes",
  "Unio crassus",
  "Unio elongatulus",
  "Unio pictorum",
  "Unio bruguierianus",
  "Unio desectus",
  "Unio ionicus",
  "Unio vicarius",
  "Vertigo angustior",
  "Vertigo moulinsiana",
  "Zerynthia polyxena"
)

######################## Species taxonomy ########################
necca_redlist_points <- st_read("../data/necca_redlist/points_invertebrates.gpkg")

## Austropotamobious torrentium is wrong! Austropotamobius torrentium is the correct
necca_redlist_species <- unique(necca_redlist_points$sci_name)
necca_redlist_gnr <- gna_verifier(necca_redlist_species,
                                     all_matches=T,
                                     data_sources=c(1,9,11,12)) #Catalogue of Life, Worms, Gbif,EOL
write_delim(necca_redlist_gnr,"../results/necca_redlist_gnr.tsv",delim="\t")


### gnr verifier to resolve species names
gnr_species_verifier <- gna_verifier(species_names_combined,
                                     all_matches=T,
                                     data_sources=c(1,9,11,12)) #Catalogue of Life, Worms, Gbif,EOL

write_delim(gnr_species_verifier, "../results/gnr_species_verifier.tsv",delim="\t")

species_gbif_taxonomy_search <- name_backbone_checklist(species_names_combined,verbose=TRUE)
write_delim(species_gbif_taxonomy_search, "../results/species_gbif_taxonomy_search.tsv",delim="\t")

species_gbif_taxonomy <- species_gbif_taxonomy_search |>
    filter(is_alternative=="FALSE")

write_delim(species_gbif_taxonomy, "../results/species_gbif_taxonomy_f.tsv",delim="\t")

#species_names <- as.character(invertebrates_all_species$SPECIES_NAME)
#gbif_species <- get_gbifid(species_names,ask=F)


######################## GBIF ########################
# GBIF retrieve data for all arthropod species that have been assessed in IUCN
### NOT run takes time. 
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
#
########################## species occurrences ##############################
# GBIF occurrences
## need to set GBIF credential in the .Renviron file
gbif_taxon_df <- gnr_species_verifier |>
    filter(dataSourceId==11) 

gbif_taxon_keys <- as.numeric(gbif_taxon_df$recordId)
gbif_taxon_ids_interactive <- get_gbifid(species_names_combined)

#gbif_taxon_keys <- as.numeric(na.omit(species_gbif_df$gbifid.ids))
#
## run once to request the download from the server
occ_download(
pred_in("taxonKey", gbif_taxon_keys), # important to use pred_in
pred("country", "GR"),
pred("hasCoordinate", TRUE),
pred("hasGeospatialIssue", FALSE),
format = "SIMPLE_CSV"
)

# download on 16/7/2025
occ_download_wait('0024201-250711103210423')
#  Status: SUCCEEDED
#  DOI: 10.15468/dl.dyfhdg
#  Format: SIMPLE_CSV
#  Download key: 0024201-250711103210423
#  Created: 2025-07-16T14:58:35.297+00:00
#  Modified: 2025-07-16T15:00:38.467+00:00
#  Download link: https://api.gbif.org/v1/occurrence/download/request/0024201-250711103210423.zip
#  Total records: 2498



# download on 30/5/2025
#occ_download_wait('0008134-250525065834625')
# DOI: 10.15468/dl.d2v94q 
# OLD download 
# GBIF.org (23 April 2025) GBIF Occurrence Download  https://doi.org/10.15468/dl.v7ruda
# Another key for the gbif download of 268 invertegrate species 
# is 0018673-241107131044228 

#gbif_species_occ <- occ_download_get('0024201-250711103210423') |>
#    occ_download_import()

#write_delim(gbif_species_occ, "../data/gbif_invertebrate_species_occ.tsv", delim="\t")

gbif_species_occ <- read_delim("../data/gbif_invertebrate_species_occ.tsv", delim="\t")

# check which is na in species. 
gbif_species_occ |> filter(is.na(species))

# 
gbif_species_occ$species <- ifelse(!is.na(gbif_species_occ$verbatimScientificName) &
                                    gbif_species_occ$verbatimScientificName=="Panaxia quadripunctaria",
                                     "Euplagia quadripunctaria",
                                     as.character(gbif_species_occ$species))

gbif_species_occ_sf <- gbif_species_occ |> 
    st_as_sf(coords=c("decimalLatitude","decimalLongitude"),
             remove=F,
             crs="WGS84") |>
    st_transform(4326) |> 
    mutate(species=if_else(verbatimScientificName=="Panaxia quadripunctaria","Euplagia quadripunctaria",species))

greece_regions_bbox <- st_as_sf(st_as_sfc(st_bbox(greece_regions)))

gbif_species_map <- ggplot() +
    geom_sf(greece_regions, mapping=aes()) +
    geom_point(gbif_species_occ_sf,
            mapping=aes(x=decimalLongitude,
                        y=decimalLatitude,
                        color=species),
            size=1.8,
            alpha=0.8,
            show.legend=T) +
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../figures/gbif_species_map.png", 
       plot=gbif_species_map, 
       height = 40, 
       width = 60,
       dpi = 300, 
       units="cm",
       device="png")

## Greece only the search in country of GBIF resolves this.

# Define the bounding box coordinates
#xmin <- 19.37359
#ymin <- 34.80202
#xmax <- 29.64306
#ymax <- 41.7485
#
#gbif_species_occ_gr <- gbif_species_occ_sf |>
#    filter(decimalLongitude > 19.37359 & decimalLongitude<29.64306,
#           decimalLatitude>34.80202 & decimalLatitude < 41.7485 )
#
#gbif_species_gr_map <- ggplot() +
#    geom_sf(greece_regions, mapping=aes()) +
#    geom_point(gbif_species_occ_gr,
#            mapping=aes(x=decimalLongitude,
#                        y=decimalLatitude,
#                        color=species),
#            size=1.8,
#            alpha=0.8,
#            show.legend=T) +
#    coord_sf(crs="WGS84") +
#    theme_bw()+
#    theme(axis.title=element_blank(),
#          axis.text=element_text(colour="black"),
#          legend.title = element_text(size=8),
#          legend.position = "bottom",
#          legend.box.background = element_blank())
#
#ggsave("../figures/gbif_species_gr_map.png", 
#       plot=gbif_species_gr_map, 
#       height = 40, 
#       width = 40,
#       dpi = 300, 
#       units="cm",
#       device="png")

################################## world Clim ####################################


####### World clim #########
world_clim_directory <- "/Users/talos/Documents/spatial_data/world_clim/wc2.1_30s_bio/"
output_directory <- "/Users/talos/Documents/programming_projects/necca_epopteia/spatial_data/world_clim_greece/"

world_clim_files <- list.files(world_clim_directory, pattern = "\\.tif$", full.names = TRUE)

for (f in world_clim_files) {
    
    if (grepl("*.tif$", f)) {
        
        #read_raster
        path_raster <- f
        crop_raster_by_shp(path_raster,greece_regions_bbox_sf,output_directory)
        
        rm(path_raster,raster_tmp,output_raster)

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
        crop_raster_by_shp(path_raster,greece_regions_bbox_sf,output_directory)
        
        rm(path_raster,raster_tmp,output_raster)

    }else{
        
        print(f, " not a tif")
        next
    }
}

### to delete
gr_1km <- sf::st_read("../spatial_data/eea_reference_grid/gr_1km.shp") 

#greece_bbox_LAEA <- st_bbox(gr_1km) |> st_as_sfc()

### corine land cover raster 100m
corine_file <- "/Users/talos/Documents/spatial_data/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif"
corine_gr_d <- "../spatial_data/u2018_clc2018_v2020_20u1_raster100m_gr/"

corine_europe <- rast(corine_file)
# crop
#crop_raster_by_shp(corine_file,greece_regions_bbox_LAEA,corine_gr_d)

corine_gr <- crop(corine_europe,greece_regions_bbox_LAEA)

#gr_1km_terra <- mask(gr_1km_c, greece_regions_bbox_sf)

writeRaster(corine_gr, "../spatial_data/u2018_clc2018_v2020_20u1_raster100m_gr/crop_U2018_CLC2018_V2020_20u1.tif",overwrite=TRUE)

corine_gr <- rast("../spatial_data/u2018_clc2018_v2020_20u1_raster100m_gr/crop_U2018_CLC2018_V2020_20u1.tif")

#### EU DEM
###https://ec.europa.eu/eurostat/web/gisco/geodata/digital-elevation-model/eu-dem#DD
# crop the eu dem file, it is 20gb.
eu_dem_file <- "/Users/talos/Documents/spatial_data/EU_DEM_mosaic_1000K/eudem_dem_3035_europe.tif"
eurodem_gr_d <- "../spatial_data/EU_DEM_mosaic_5deg_gr/"

eu_dem_e <- rast(eu_dem_file)

eu_dem_gr <- crop(eu_dem_e,greece_regions_bbox_LAEA)

terra::writeRaster(eu_dem_gr,
                   "../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_3035_europe.tif",
                   overwrite=TRUE)

#crop_raster_by_shp(eu_dem_file,greece_regions_bbox_sf,eurodem_gr_d)


# Then reload in different name and change crs
eu_dem <- rast("../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_4258_gr.tif")

eu_dem_wgs84 <- terra::project(eu_dem, "EPSG:4326")


#### EU dem slope
eu_dem_slope_file <- "/Users/talos/Documents/spatial_data/EUD_CP_SLOP_mosaic/eudem_slop_3035_europe.tif"

eu_dem_slope_gr_d <- "../spatial_data/EU_DEM_slope_gr/"

eu_dem_slope_e <- rast(eu_dem_slope_file)

eu_slope_gr <- crop(eu_dem_slope_e,greece_regions_bbox_LAEA)

terra::writeRaster(eu_slope_gr,
                   "../spatial_data/EU_DEM_slope_gr/crop_eudem_slop_3035_europe.tif",
                   overwrite=TRUE)
#crop_raster_by_shp(eu_dem_slope_file,greece_regions_bbox_sf,eu_dem_slope_gr_d)

# Then reload in different name and change crs
# the conversion happened succesfully in the function
#eu_dem_slope_wgs84 <- terra::project(eu_dem_slope, "EPSG:4326")

#terra::writeRaster(eu_dem_slope_wgs84,
#                   "../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_slope_4326_gr.tif",
#                   overwrite=TRUE)

### Ecosystem types of Europe 3.1
ecosystem_types <- "/Users/talos/Documents/spatial_data/Ecosystem_types_of_Europe_version_3.1_Terrestrial/eea_r_3035_100_m_etm-terrestrial-c_2012_v3-1_r00.tif"

ecosystem_types_gr_d <- "../spatial_data/ecosystem_types_gr/"
crop_raster_by_shp(ecosystem_types,greece_regions_bbox_sf,ecosystem_types_gr_d)


#ecosystem_types_gr <- rast("/Users/talos/Documents/programming_projects/necca_epopteia/spatial_data/ecosystem_types_gr/crop_eea_r_3035_100_m_etm-terrestrial-c_2012_v3-1_r00.tif")

############### rasterize ###############
### reference grid ###
#load
gr_1km <- sf::st_read("../spatial_data/eea_reference_grid/gr_1km.shp") 

grid_sf <- gr_1km

# 2) Make/ensure a unique cell ID (character or integer)
#if (!"CELLCODE" %in% names(grid_sf)) {
#  grid_sf$cell_id <- seq_len(nrow(grid_sf))
#  grid_v <- vect(grid_sf)
#}

# 3) Build a raster template aligned to the grid
#    Set resolution to 1000 for 1-km grid (or 10000 for 10-km)

r_tmpl <- rast(grid_v, res = 1000)         # use 10000 for 10 km

origin(r_tmpl) <- c(0, 0)                  # snap to EEA grid origin
values(r_tmpl) <- NA
crs(r_tmpl) <- crs(grid_v)

# 4) Rasterize the grid polygons to cell IDs (categorical raster)
r_id <- rasterize(grid_v, r_tmpl, field = "CELLCODE")


# Save
writeRaster(r_tmpl, "../spatial_data/eea_reference_grid/EEA_grid_id_1km.tif", overwrite = TRUE, datatype = "INT4S")

# rasterize
polygon_sf <- gr_1km 

raster_template <- rast(ext = ext(polygon_sf), resolution = 1000)  # 1000 meters = 1 km resolution

raster_object <- rasterize(polygon_sf, raster_template, field = "CELLCODE", background = NA)

gr_1km_rast <- raster_object
crs(gr_1km_rast) <- "EPSG:3035"
writeRaster(raster_object, "../spatial_data/eea_reference_grid/gr_1km.tif",overwrite=TRUE)

# Reproject using resample to align with the polygon's CRS

gr_1km_rast_wgs <- project(gr_1km_rast,"EPSG:4326")

writeRaster(gr_1km_rast_wgs, "../spatial_data/eea_reference_grid/gr_1km_wgs.tif",overwrite=TRUE)

### keep only land

bbox_sf <- st_as_sf(st_as_sfc(st_bbox(greece_regions)))
gr_1km_c <- crop(gr_1km_rast_wgs,greece_regions_bbox_sf)
gr_1km_terra <- mask(gr_1km_c, greece_regions_bbox_sf)

writeRaster(gr_1km_terra, "../spatial_data/eea_reference_grid/gr_1km_terra.tif",overwrite=TRUE)

### vegetation map


vegetation_map <- sf::st_read("../spatial_data/Vegetation_map_Greece/D_xabxg_VPG_60-98_GEO.shp",
                              options = "ENCODING=WINDOWS-1253")

vegetation_map_wgs <- st_transform(vegetation_map,4326) |>
    st_make_valid()

veg_data <- data.frame(
  A_VEG_TYPE = c("ΕΛΑ", "ΕΡΛ", "ΠΜΑ", "ΠΛΔ", "ΠΔΑ", "ΠΧΑ", "ΠΚΟ", "ΠΘΑ", "ΚΠΡ",
                 "ΑΡΚ", "ΟΞΥ", "ΔΡΥ", "ΚΑΣ", "ΣΗΜ", "ΣΦΕ", "ΦΙΛ", "ΦΠΛ", "ΠΑΡ",
                 "ΕΥΚ", "ΦΟΙ", "ΘΑΜ", "ΦΘΑ", "ΛΙΒ", "ΑΓΟ", "ΟΙΚ", "ΓΚΑ", "ΓΚΕ",
                 "ΛΧΡ", "ΛΙΜ"),
  A_VEG_NAME = c("Ελάτη", "Ερυθρελάτη", "Πεύκη μαύρη", "Πεύκη λευκόδερμη", "Πεύκη δασική",
                 "Πεύκη χαλέπιος", "Πεύκη κουκουναριά", "Πεύκη θαλασσία", "Κυπαρίσσι",
                 "Άρκευθος", "Οξυά", "Δρύς", "Καστανιά", "Σημύδα", "Σφένδαμος",
                 "Φιλύρα", "Φυλλοβόλα πλατύφυλλα", "Παραποτάμια βλάστηση", "Ευκάλυπτος",
                 "Φοίνικας", "Θάμνοι", "Φυλλοβόλοι θάμνοι", "Λιβάδια, αραιά ξυλ. βλάστηση",
                 "Άγονα", "Οικισμοί", "Γεωργ. καλλιέργειες", "Γεωργ. καλλιέργειες εγκατ.",
                 "Λοιπές χρήσεις", "Λίμνη"),
  stringsAsFactors = FALSE
)

veg_rast <- vect(vegetation_map)

veg_rast <- project(veg_rast, crs(slope_rast))

veg_rasterized <- rasterize(veg_rast, slope_rast, field = "A_VEG_TYPE")


writeRaster(veg_rasterized, "../spatial_data/Vegetation_map_Greece/D_xabxg_VPG_60-98_GEO.tif")
