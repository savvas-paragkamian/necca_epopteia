#!/usr/bin/env Rscript

## Script name: spatial_analysis.R
##
## Purpose of script:
## 1. Load all data
## 2. Enrich occurrences with spatial info
## 3. Quality control based on custom filters (e.g elevation) of species
## 4. Compose a sample presence file file with all relavant information
## 5. Assign filters in includeDistribution using
##      a) data origin (datasetName)
##      b) spatial info (distance from borders and elevation)
##      c) species specific curation
## 6. Calculations
##      a) species metrics, population, population_n2k, distribution, range
## 7. Maps
##
## Author: Savvas Paragkamian
## 
## Runtime: 2 minutes without creating maps
##
## Date Created: 2024-11-06

library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(units)
library(ggpubr)
library(ggnewscale)
source("necca_spatial_functions.R")

#----------------------------------------------------------------------------#
########################### Load Spatial Data ################################
#----------------------------------------------------------------------------#

### greece 
greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")

greece_regions_ETRS89 <- st_transform(greece_regions, 3035)
### EEA reference grid
eea_1km <- sf::st_read("../spatial_data/eea_reference_grid/gr_1km.shp")
eea_1km_wgs <- st_transform(eea_1km,4326)
eea_1km_ETRS89 <- st_transform(eea_1km, 3035)

eea_1km_ETRS89_min <- eea_1km_ETRS89 |>
    dplyr::select(CELLCODE,geometry) |>
    rename("CELLCODE_eea_1km"="CELLCODE")

eea_10km <- sf::st_read("../spatial_data/eea_reference_grid/gr_10km.shp")
eea_10km_wgs <- st_transform(eea_10km,4326)
eea_10km_ETRS89 <- st_transform(eea_10km, 3035)

eea_10km_ETRS89_c <- st_centroid(eea_10km_ETRS89)
eea_10km_ETRS89_c_wgs <- st_transform(eea_10km_ETRS89_c,4326)
### eea 10km only land and in greece
grids_in_greece <- st_intersection(eea_10km_ETRS89, greece_regions_ETRS89)

# remove not used columns and calculate area
grids_in_greece_min <- grids_in_greece |>
    dplyr::select(CELLCODE,geometry) |>
    rename("CELLCODE_eea_10km"="CELLCODE") |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    group_by(CELLCODE_eea_10km) |>
    summarise(AreaKm2=sum(area))

### Natura2000

N2000_v32 <- sf::st_read("../spatial_data/N2000_spatial_GR_2021_12_09_v32/N2000_spatial_GR_2021_12_09_v32.shp")

#transform Natura2000
N2000_v32_ETRS89 <- st_transform(N2000_v32, 3035) |>
    dplyr::select(SITECODE,SITETYPE,geometry)

#wgs
#N2000_v32_wgs <- st_transform(N2000_v32,4326)
#N2000_v32_wgs_sci <- N2000_v32_wgs |>
#    filter(SITETYPE=="SCI")

#N2000_v32_wgs_spa <- N2000_v32_wgs |>
#    filter(SITETYPE=="SPA")

#N2000_v32_wgs_scispa <- N2000_v32_wgs |>
#    filter(SITETYPE=="SCISPA")

#
#----------------------------------------------------------------------------#
########################### Load Species Point Data ##########################
########################### derived from species_enrichment.R ################
#----------------------------------------------------------------------------#
### Taxonomy
species_taxonomy <- read_delim("../results/species_gbif_taxonomy_curated.tsv",delim="\t")

species_taxonomy_only <- species_taxonomy |>
    dplyr::select(species,verbatim_name)

species_names_combined <- as.character(species_taxonomy$verbatim_name)

### Species occurrences integrated ###
species_samples_art17_all <- read_delim("../results/species_samples_art17_all.tsv", delim="\t")
# species_occurrences_invertebrates <- read_delim("../results/species_occurrences_invertebrates.tsv",delim="\t")

## remove points without coordinates
### distinct all points

species_samples_simple <- species_samples_art17_all |> 
    filter(!is.na(decimalLatitude)) |>
    distinct(submittedName,
             decimalLatitude,
             decimalLongitude,
             datasetName,
             basisOfRecord,
             individualCount,
             species) |>
    #left_join(species_taxonomy_only, by=c("submittedName"="verbatim_name")) |>
    mutate(species = gsub("Hirudo medicinalis","Hirudo verbana",species)) |>
    mutate(species = gsub("Unio vicarius","Unio crassus",species)) |>
    mutate(species = gsub("Unio ionicus","Unio crassus",species)) |>
    mutate(species = gsub("Unio bruguierianus","Unio crassus",species)) |>
    mutate(species = gsub("Unio desectus","Unio crassus",species)) |>
    mutate(species = gsub("Polyommatus eros","Polyommatus eroides",species)) |>
    mutate(species = gsub("Paracossulus thrips","Catopta thrips",species)) |>
    mutate(species = gsub("Hirudo medicinalis","Hirudo verbana",species)) |> 
    mutate(species = gsub("Osmoderma eremita","Osmoderma lassallei",species))

# shapefile for visualisation
points_sf <- species_samples_simple |> 
    filter(!is.na(decimalLatitude)) |>
    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),
             remove=F,
             crs="WGS84")

    # find the duplicates
points_duplicates <- species_samples_art17_all |>
    filter(!is.na(decimalLatitude)) |>
    group_by(across(everything())) |>
    summarise(n = n(), .groups = "drop") |>
    filter(n>1)

## sanity check
print("is the number of rows the same with duplicates and removal of NAs in coordinates?")
nrow(species_samples_art17_all |> filter(!is.na(decimalLatitude)) ) == nrow(points_sf) + sum(points_duplicates$n) - nrow(points_duplicates)

# tranform to ETRS89
points_sf_ETRS89 <- st_transform(points_sf, 3035)

# add the cells 
points_sf_cells <- st_join(points_sf_ETRS89,eea_10km_ETRS89)

#----------------------------------------------------------------------------#
## point samples is ready to use for the rest of the analysis
species_samples <- points_sf_cells |>
    st_drop_geometry() |>
    rename("CELLCODE_eea_10km"="CELLCODE") |>
    distinct(submittedName,
             decimalLatitude,
             decimalLongitude,
             datasetName,
             basisOfRecord,
             individualCount,
             CELLCODE_eea_10km,
             species)

#----------------------------------------------------------------------------#
#################### Distribution data from previous report ##################
######################### National report_2013_2018_shp ######################
#----------------------------------------------------------------------------#
### this dataset contains multipolygons, one for each species.
### each polygon from this multipolygon is a grid from the EEA 10X10 reference
### grid. The cellcode is missing!
species_dist_national_rep <- sf::st_read("../spatial_data/National report_2013_2018_shp/GR_Art17_species_distribution.shp")
species_dist_national_rep_sens <- sf::st_read("../spatial_data/National report_2013_2018_shp/GR_Art17_species_distribution_sensitive.shp")

## merge the two distributions

# change the names of the species
# one feature in invalid
# st_is_valid(species_dist_national_rep)
# st_is_valid(species_dist_national_rep_sens)
species_dist_national_E1 <- rbind(species_dist_national_rep,species_dist_national_rep_sens) |>
    st_make_valid() |>
    mutate(submittedName=gsub("\\* ?","",iname)) |>
    mutate(submittedName=gsub(" Complex","",submittedName)) |>
    mutate(submittedName = gsub("Osmoderma eremita","Osmoderma lassallei",submittedName)) |>
    filter(submittedName %in% species_names_combined) |>
    left_join(species_taxonomy_only, by=c("submittedName"="verbatim_name")) |>
    mutate(species = gsub("Osmoderma eremita","Osmoderma lassallei",species)) |>
    mutate(species = gsub("Hirudo medicinalis","Hirudo verbana",species)) |>
    mutate(species = gsub("Unio vicarius","Unio crassus",species)) |>
    mutate(species = gsub("Unio ionicus","Unio crassus",species)) |>
    mutate(species = gsub("Unio bruguierianus","Unio crassus",species)) |>
    mutate(species = gsub("Unio desectus","Unio crassus",species)) |>
    mutate(species = gsub("Polyommatus eros","Polyommatus eroides",species)) |>
    mutate(species = gsub("Paracossulus thrips","Catopta thrips",species)) 

## add the Cellcode with st_within join to have only those that are
## from the same grids
species_dist_national_E1_eea <- st_join(eea_10km_ETRS89,
                                        species_dist_national_E1,
                                        join = st_within, 
                                        left = FALSE,
                                        largest = FALSE)

st_write(species_dist_national_E1_eea, "../results/species_dist_national_E1_eea/species_dist_national_E1_eea.shp", layer = "species_dist_national_E1_eea", delete_layer = TRUE)


st_write(species_dist_national_E1_eea, "../results/species_dist_national_E1_eea.gpkg", layer = "species_dist_national_E1_eea", delete_layer = TRUE)

# centroids
species_dist_national_E1_centroids <- st_centroid(species_dist_national_E1_eea)
species_dist_national_E1_coords <- st_coordinates(species_dist_national_E1_centroids)

# Combine metadata + coordinates into a data.frame
species_dist_national_E1_eea_coords <- cbind(
                                             species_dist_national_E1_eea,
                                             decimalLongitude = species_dist_national_E1_coords[, 1],
                                             decimalLatitude = species_dist_national_E1_coords[, 2])


# transform the columns
species_dist_national_E1_eea <- species_dist_national_E1_eea_coords |>
    distinct(submittedName, species,CELLCODE,decimalLongitude,decimalLatitude) |>
    rename("CELLCODE_eea_10km"="CELLCODE") |>
    mutate(
           datasetName="DistrMap_2013_2018",
           basisOfRecord="ESTIMATED_CENTROID",
           DistrMap_2013_2018=TRUE
           ) 

species_dist_national_E1_eea_coords_laea <- species_dist_national_E1_eea |>
    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),
             remove=T,
             crs=crs(species_dist_national_rep))

species_dist_national_E1_eea_coords_wgs <-  st_transform(species_dist_national_E1_eea_coords_laea,4326)

species_dist_national_E1_coords_w_c <- st_coordinates(species_dist_national_E1_eea_coords_wgs)

# Combine metadata + coordinates into a data.frame
species_dist_national_E1_eea_coords_wgs_f <- cbind(
                                             species_dist_national_E1_eea_coords_wgs,
                                             decimalLongitude = species_dist_national_E1_coords_w_c[, 1],
                                             decimalLatitude = species_dist_national_E1_coords_w_c[, 2])

species_dist_national_E1_eea_coords_wgs <- species_dist_national_E1_eea_coords_wgs_f

st_write(species_dist_national_E1_eea_coords_wgs, "../results/species_dist_national_E1_eea_coords_wgs.gpkg", layer = "species_dist_national_E1_eea_coords_wgs", delete_layer = TRUE)

#----------------------------------------------------------------------------#
########## Parnassius apollo distribution data Action Plan 2019 #############
#----------------------------------------------------------------------------#

# P apollo
#----------------------------------------------------------------------------#
p_apollo_dist <- st_read("../data/Parnassius apollo AP 2019/AP_Papollo_Distribution_LAEA.shp") |>
    mutate(
           datasetName="Action Plan 2019",
           basisOfRecord="ESTIMATED_CENTROID"
           ) 
# centroids
p_apollo_centroids <- st_centroid(p_apollo_dist)
p_apollo_coords <- st_coordinates(p_apollo_centroids)

# Combine metadata + coordinates into a data.frame
p_apollo_eea_coords <- cbind(
                                             p_apollo_dist,
                                             decimalLongitude = p_apollo_coords[, 1],
                                             decimalLatitude = p_apollo_coords[, 2])



p_apollo_coords_laea <- p_apollo_eea_coords |>
    st_drop_geometry() |>
    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),
             remove=T,
             crs=crs(eea_10km_ETRS89))

p_apollo_coords_wgs <-  st_transform(p_apollo_coords_laea,4326)

p_apollo_coords_d <- st_coordinates(p_apollo_coords_wgs)

p_apollo_eea_coords_wgs <- cbind(p_apollo_coords_wgs,
                                             decimalLongitude = p_apollo_coords_d[, 1],
                                             decimalLatitude = p_apollo_coords_d[, 2])

## prepare the df of P apollo from Action Plan 2019
p_apollo_dist_d <- p_apollo_eea_coords_wgs |>
    st_drop_geometry() |>
    dplyr::select(-c(EOFORIGIN,NOFORIGIN)) |>
    distinct() |>
    rename("CELLCODE_eea_10km"="CELLCODE") |>
    mutate(
           species="Parnassius apollo",
           individualCount=NA,
           submittedName="Parnassius apollo"
           ) 

#----------------------------------------------------------------------------#
############################# Presence File Compilation #####################
######### Merging previous files
#----------------------------------------------------------------------------#

# input is species_samples
# add dist of p apollo
species_samples_eea_minimum <- species_samples |>
    rbind(p_apollo_dist_d)

# national distribution as dataset
species_dist_national_dataset <- species_dist_national_E1_eea_coords_wgs |>
    st_drop_geometry() |>
    dplyr::select(-c(DistrMap_2013_2018))

# find which species-cellcode pairs are also in the report DistrMap_2013_2018
#
species_dist_national_minimum <- species_dist_national_E1_eea_coords_wgs |>
    distinct(species,CELLCODE_eea_10km,DistrMap_2013_2018)

# merge the column of DistrMap_2013_2018
species_samples_eea_minimum_dist <- species_samples_eea_minimum |>
    left_join(species_dist_national_minimum,
              by=c("species"="species","CELLCODE_eea_10km"="CELLCODE_eea_10km")) |>
    mutate(DistrMap_2013_2018 = if_else(is.na(DistrMap_2013_2018), FALSE, DistrMap_2013_2018)) |>
    mutate(composite_key = paste0(species,"_",CELLCODE_eea_10km))

# remove gbif from the calculation of orphans cells
species_samples_eea_minimum_dist_no_gbif <- species_samples_eea_minimum_dist |>
    filter(datasetName!="GBIF")

############### identify orphans cells i.e cells without point data##########
# merge the national report as a dataset
species_dist_national_comp <- species_dist_national_dataset |>
    mutate(composite_key = paste0(species,"_",CELLCODE_eea_10km))

national_orphans <- species_dist_national_comp |>
    filter(!c(composite_key %in% unique(species_samples_eea_minimum_dist_no_gbif$composite_key))) |>
    mutate(DistrMap_2013_2018 = TRUE) |>
    mutate(individualCount = NA) #|>
    #mutate(includeDistribution = TRUE) 

##### evaluate the origin of orphans from the original file

E1X_DB_ref_samples_data <- read_xlsx("../data/Ε1Χ_ΒΔ_ΒΙΒΛΙΟΓΡΑΦΙΑΣ_ΑΣΠ_20250802.xlsx",
                                    sheet="Εξάπλωση ειδών και τ.ο.",
                                    col_names=T) |> slice(-1) |> filter(!is.na(`Κωδικός Αναφοράς`))

E1X_DB_refs_data <- read_xlsx("../data/Ε1Χ_ΒΔ_ΒΙΒΛΙΟΓΡΑΦΙΑΣ_ΑΣΠ_20250802.xlsx",
                                    sheet="Βιβλιογραφία",
                                    col_names=T) |> slice(-1) |> filter(!is.na(`Κωδικός Αναφοράς`))


E1X_DB_ref_samples_data$decimalLatitude <- as.numeric(E1X_DB_ref_samples_data$`Γεωγραφικό Πλάτος (WGS84)`)
E1X_DB_ref_samples_data$decimalLongitude <- as.numeric(E1X_DB_ref_samples_data$`Γεωγραφικό Μήκος (WGS84)`)

E1X_DB_ref_all <- E1X_DB_ref_samples_data |>
    left_join(E1X_DB_refs_data,
              by=c("Κωδικός Αναφοράς"="Κωδικός Αναφοράς")) |>
    mutate(datasetName="E1X_DB_references") |>
    mutate(basisOfRecord="MaterialCitation") |>
    mutate(submittedName=`Ονομασία είδους`) |>
    mutate(individualCount=as.numeric(`Πλήθος ατόμων`)) |>
    mutate(CELLCODE_eea_10km = `Κωδικός Κελιού Πλέγματος 10x10km (ETRS LAEA)`) |>
    dplyr::select(datasetName,
                  submittedName,
                  individualCount,
                  CELLCODE_eea_10km,
                  decimalLatitude,
                  decimalLongitude
                  ) |>
    left_join(species_taxonomy, by=c("submittedName"="verbatim_name"))

E1X_DB_ref_all_cells <- E1X_DB_ref_all |>
    filter(is.na(decimalLatitude)) |>
    distinct(datasetName,individualCount,CELLCODE_eea_10km,species) |>
    mutate(composite_key = paste0(species,"_",CELLCODE_eea_10km))

E1X_DB_ref_all_cells_orphans <- E1X_DB_ref_all_cells[which(E1X_DB_ref_all_cells$composite_key %in% national_orphans$composite_key),]

## orphans that 
orphans_in_E1X <- national_orphans[which(!(national_orphans$composite_key %in% unique(E1X_DB_ref_all_cells$composite_key))),] |> as_tibble() |> distinct(composite_key)

## the reverse, E1X in orphan cells
E1X_in_orphans <- E1X_DB_ref_all_cells[which(!(E1X_DB_ref_all_cells$composite_key %in% unique(national_orphans$composite_key))),] |> as_tibble() |> distinct(composite_key)

#----------------------------------------------------------------------------#
# final file with the previous distribution and orphan cells
# compose final presence file
# some species filtering are applied to inform the calculation of Distribution
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
## presence file with orphans
species_samples_presence_minimum <- species_samples_eea_minimum_dist |>
    rbind(national_orphans) |>
    dplyr::select(-composite_key) 

#----------------------------------------------------------------------------#
###################### Extract data to species occurrences ####################
#----------------------------------------------------------------------------#

################ elevation ################

### EU DEM Greece
eu_dem_gr <- rast("../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_4258_europe.tif")
#eu_dem_gr_4258 <- rast("../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_4258_gr.tif")

species_samples_presence_minimum_sf <- species_samples_presence_minimum |>
    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),
             remove=F,
             crs="WGS84")

vect_points <- terra::vect(species_samples_presence_minimum_sf)
  
vals <- terra::extract(eu_dem_gr, vect_points)[, -1, drop = FALSE]

species_samples_presence_sf <- cbind(species_samples_presence_minimum_sf, vals)

no_elevation <- species_samples_presence_sf[which(is.na(species_samples_presence_sf$eudem_dem_4258_europe)),]

print("how many points don't have elevation?")
nrow(no_elevation)
#################### distance from land and/or border of Greece #############

species_samples_presence_sf_ETRS89 <- st_transform(species_samples_presence_sf, 3035)

dist_matrix <- st_distance(species_samples_presence_sf_ETRS89, greece_regions_ETRS89)
species_samples_presence_sf_ETRS89$minimumDistanceFromBorders <- apply(dist_matrix, 1, min)

################### analysis of points distance from land ##############
points <- species_samples_presence_sf_ETRS89

# points inside
points_inside_or_touching <- points[points$minimumDistanceFromBorders ==0, ]

# Keep only points within 500, 2000 m of any polygon
points_away_500m <- points[points$minimumDistanceFromBorders >500, ]
# Keep only points with distance > 2000 of any polygon
points_away_2km <- points[points$minimumDistanceFromBorders > 2000, ]
# Keep all points within 500 m of any polygon
points_500m <- points[points$minimumDistanceFromBorders > 0 & points$minimumDistanceFromBorders <500, ]

## save points outside
points_outside_cols <- points |>
    filter(minimumDistanceFromBorders >0) |>
    dplyr::select(species, datasetName,basisOfRecord,geometry,minimumDistanceFromBorders,eudem_dem_4258_europe) |>
    arrange(species)

write_delim(points_outside_cols,
            "../results/species_occurrences_outside_gr_land.tsv",
            delim="\t")
# plot

points_gr_plot <- ggplot() +
    geom_sf(data = eea_10km_ETRS89, aes(color = "eea"), fill="transparent", shape = 16, show.legend = TRUE) +
    geom_sf(data = greece_regions_ETRS89, aes(color = "Polygons"), fill = NA, show.legend = TRUE) +
    geom_sf(data = points_inside_or_touching, aes(color = "Points inside"), shape = 1, show.legend = TRUE) +
    geom_sf(data = points_500m, aes(color = "Points away < 500m"), shape = 1, show.legend = TRUE) +
    geom_sf(data = points_away_500m, aes(color = "Points away > 500m"), shape = 1, show.legend = TRUE) +
    scale_color_manual(
      name = "Feature Type",
      values = c(
        "Gr borders" = "blue",
        "Points inside" = "black",
        "Points away > 500m" = "blue",
        "Points away < 500m" = "red"
      )
    ) +
    theme_bw()+
    theme(legend.position="inside",
          legend.position.inside = c(0.87, 0.75))

ggsave(plot=points_gr_plot,
       "../figures/map_occurrences_borders_filtering.png",
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")

#----------------------------------------------------------------------------#
####################### Custom filters for includeDistribution################
##############################################################################
#----------------------------------------------------------------------------#

species_samples_presence_sf_ETRS89_dist <- species_samples_presence_sf_ETRS89 |>
    mutate(species = if_else(species=="Zerynthia polyxena" & decimalLatitude < 36,
                                    "Zerynthia cretica",
                                    species)) |>
    mutate(includeDistribution=TRUE) |>
    mutate(includeDistribution=if_else(datasetName=="GBIF",FALSE,includeDistribution) ) |>
    mutate(includeDistribution=if_else(datasetName=="E1X_DB_references" & DistrMap_2013_2018==FALSE,
                                       FALSE,
                                       includeDistribution) ) |>
    mutate(includeDistribution=if_else(species=="Zerynthia cretica", 
                                       FALSE,
                                       includeDistribution)) |>
    mutate(includeDistribution=if_else(species=="Vertigo moulinsiana",
                                       FALSE,
                                       includeDistribution)) |>
    mutate(includeDistribution=if_else(minimumDistanceFromBorders > 500 &
                                       basisOfRecord!="ESTIMATED_CENTROID",
                                       FALSE,
                                       includeDistribution)) |>
    mutate(includeDistribution=if_else(species=="Stenobothrus eurasius" & CELLCODE_eea_10km %in% c("10kmE535N184"),
                                       FALSE,
                                       includeDistribution)) |>
    mutate(includeDistribution=if_else(species=="Probaticus subrugosus" & CELLCODE_eea_10km %in% c("10kmE524N175"),
                                       FALSE,
                                       includeDistribution)) |>
    mutate(includeDistribution=if_else(species=="Pseudophilotes bavius" & CELLCODE_eea_10km %in% c("10kmE535N176"),
                                       FALSE,
                                       includeDistribution)) |>
    mutate(includeDistribution=if_else(species == "Parnassius apollo" & 
                                       eudem_dem_4258_europe <= 450 & 
                                       datasetName!="Action Plan 2019",
                                       FALSE,
                                       includeDistribution)) |>
    mutate(includeDistribution=if_else(species=="Parnassius apollo" & datasetName %in% c("E1X_DB_references","DistrMap_2013_2018"),
                                       FALSE,
                                       includeDistribution)) |>
#    mutate(includeDistribution=if_else(species=="Parnassius apollo" & 
#                                       minimumDistanceFromBorders > 100 & ## remove after the correct altitute file inclusion
#                                       datasetName=="E1X_MDPP_2014_2024",
#                                   TRUE,
#                                   includeDistribution)) |>
    mutate(includeDistribution=if_else(species=="Rosalia alpina" & 
                                       minimumDistanceFromBorders > 500 &
                                       datasetName=="NECCA_redlist",
                                   TRUE,
                                   includeDistribution))


# includePopulation filtering
species_samples_presence_eea_1km <- st_join(species_samples_presence_sf_ETRS89_dist,
                      eea_1km_ETRS89_min)

species_samples_presence_pop <- species_samples_presence_eea_1km |>
    mutate(includePopulation=if_else(datasetName %in% c(
                                                        "E1X_DB_references",
                                                        "DistrMap_2013_2018"),
                                     FALSE,
                                     includeDistribution)) |>
    mutate(includePopulation=if_else(species=="Vertigo angustior",
                                     FALSE,
                                     includePopulation)) 

#----------------------------------------------------------------------------#
########## FINAL FILE WITH ALL PRESENCE INFORMATION PER SPECIES #############
#----------------------------------------------------------------------------#

species_samples_presence_final_sf <- species_samples_presence_pop |>
    filter(datasetName!="Invertebrates_records_Olga") 

species_samples_presence_final <- species_samples_presence_pop |>
    st_drop_geometry() |>
    filter(datasetName!="Invertebrates_records_Olga") 

write_delim(species_samples_presence_final,
            "../results/species_samples_presence_final.tsv",
            delim="\t")

#----------------------------------------------------------------------------#
#### Distribution ####
#----------------------------------------------------------------------------#

distributions <- species_samples_presence_final |>
    left_join(grids_in_greece_min, by=c("CELLCODE_eea_10km"="CELLCODE_eea_10km")) |>
    filter(includeDistribution==TRUE) |>
    distinct(species, CELLCODE_eea_10km,AreaKm2) |>
    group_by(species) |>
    summarise(CellCount=n(), AreaKm2=round(sum(AreaKm2),digits=2))

write_delim(distributions, "../results/distributions_presence_final.tsv",delim="\t")

#----------------------------------------------------------------------------#
## population
#----------------------------------------------------------------------------#
populations_all <- species_samples_presence_final |>
    filter(includeDistribution==TRUE) |>
    #filter(includePopulation==TRUE) |>
    distinct(species, CELLCODE_eea_1km) |>
    group_by(species) |>
    summarise(n_1km=n(),
              .groups="keep") |>
    mutate(method="with_refs") |>
    ungroup()

populations_no_ref <- species_samples_presence_final |>
    filter(includePopulation==TRUE) |>
    filter(datasetName!="E1X_DB_references") |>
    distinct(species, CELLCODE_eea_1km) |>
    group_by(species) |>
    summarise(n_1km=n(),
              .groups="keep") |>
    mutate(method="no_refs") |>
    ungroup()

populations <- populations_all |>
    rbind(populations_no_ref) |>
    pivot_wider(names_from=method,values_from=n_1km)

write_delim(populations, "../results/populations_presence_final.tsv",delim="\t")
#----------------------------------------------------------------------------#
## Natura2000
#----------------------------------------------------------------------------#

# include Natura2000 columns : SITECODE and SITETYPE
#N2000_v32_ETRS89

species_samples_presence_nat <- st_join(species_samples_presence_final_sf,
                      N2000_v32_ETRS89)

## population natura all for table_12_1
populations_nat_table_12_1 <- species_samples_presence_nat |>
    st_drop_geometry() |>
    filter(includePopulation==TRUE) |>
    #filter(includeDistribution==TRUE) |>
    filter(SITETYPE %in% c("SCI","SCISPA")) |>
    distinct(species, CELLCODE_eea_1km) |>
    group_by(species) |>
    summarise(n_1km=n(),
              .groups="keep") |>
    mutate(method="no_refs_no_orphan") |>
    ungroup()

## population natura all from includeDistribution add species annex column
populations_nat_annex <- species_samples_presence_nat |>
    st_drop_geometry() |>
    filter(includePopulation==TRUE) |>
    filter(!(is.na(SITETYPE))) |>
    distinct(species, CELLCODE_eea_1km,SITETYPE,SITECODE) |>
    group_by(species,SITETYPE, SITECODE) |>
    summarise(n_1km=n(),
              .groups="keep") |>
    mutate(method="includePopulation") |>
    ungroup()

## populations per SITECODE using includeDistribution
populations_nat_no_refs_orphans <- species_samples_presence_nat |>
    st_drop_geometry() |>
    filter(includeDistribution==TRUE) |>
    filter(!(is.na(SITETYPE))) |>
    #filter(includePopulation==TRUE) |>
    distinct(species, CELLCODE_eea_1km,SITETYPE,SITECODE) |>
    group_by(species,SITETYPE, SITECODE) |>
    summarise(n_1km=n(),
              .groups="keep") |>
    mutate(method="includeDistribution") |>
    ungroup()

## populations per SITECODE using only the unique cells that are 
## derived only from E1X_DB_references
populations_nat_refs <- species_samples_presence_nat |>
    st_drop_geometry() |>
    filter(includeDistribution==TRUE) |>
    filter(!(is.na(SITETYPE))) |>
    distinct(species, CELLCODE_eea_1km,SITETYPE,SITECODE,datasetName) |>
    group_by(species, CELLCODE_eea_1km,SITETYPE,SITECODE) |>
    mutate(n_datasets=n()) |>
    filter(datasetName=="E1X_DB_references" & n_datasets==1) |>
    distinct(species, CELLCODE_eea_1km,SITETYPE,SITECODE) |>
    group_by(species,SITETYPE, SITECODE) |>
    summarise(n_1km=n(),
              .groups="keep") |>
    mutate(method="E1X_DB_references") |>
    ungroup()

## populations per SITECODE using only the unique cells that are 
## derived only from DistrMap_2013_2018 (orphan cells)
populations_nat_orphans <- species_samples_presence_nat |>
    st_drop_geometry() |>
    filter(includeDistribution==TRUE) |>
    filter(!(is.na(SITETYPE))) |>
    filter(datasetName=="DistrMap_2013_2018") |>
    distinct(species, CELLCODE_eea_1km,SITETYPE,SITECODE) |>
    group_by(species,SITETYPE, SITECODE) |>
    summarise(n_1km=n(),
              .groups="keep") |>
    mutate(method="DistrMap_2013_2018") |>
    ungroup()

## combine the above populations per SITETYPE
populations_nat_c <- rbind(populations_nat_annex,
                           populations_nat_no_refs_orphans,
                           populations_nat_orphans,
                           populations_nat_refs)

## composite file for filling and informing the 3.2
populations_nat_c_t <- populations_nat_c |>
    pivot_wider(names_from=method,values_from=n_1km,values_fill=0) |>
    mutate(only_from_refs=if_else(includeDistribution==E1X_DB_references,TRUE,FALSE)) |>
    mutate(only_from_refs_orphan=if_else(includeDistribution==E1X_DB_references+DistrMap_2013_2018 | includeDistribution==E1X_DB_references,
                                         TRUE,
                                         FALSE))

write_delim(populations_nat_c_t,"../results/populations_nat_c_t.tsv",delim="\t")

stop("Spatial analysis and Master Files are ready. Exiting script.")

#----------------------------------------------------------------------------#
## presence file summary
#----------------------------------------------------------------------------#


#----------------------------------------------------------------------------#
################################## MAPS ########################################
################################################################################
#----------------------------------------------------------------------------#

species_with_data <- unique(species_samples_presence_final_sf$species)
datasets_colors <- c(
                     "GBIF"="seagreen",
                     "NECCA_redlist"="#B31319",
                     "E1X_MDPP_2014_2024"="#FDF79C",
                     "E1X_DB"="#2BA09F",
                     "E2X_DB"="maroon1",
                     "db_refs_additional_2025"="green1",
                     "E1X_DB_references"="#141D43",
                     "Lopes-Lima et al., 2024"="#E69F00",
                     "DistrMap_2013_2018"="darkorchid1",
                     "Action Plan 2019"="#F85C29"
                     )

# base plot with Natura2000 areas of Greece
natura_colors <- c(
                   "SCI"="#E69F00",
                   "SPA"="#56B4E9",
                   "SCISPA"="#CC79A7"
)

## natura2000
g_base_n2000 <- ggplot()+
    geom_sf(greece_regions_ETRS89, mapping=aes()) +
    geom_sf(N2000_v32_ETRS89, mapping=aes(fill=SITETYPE),
            alpha=0.3,
            #colour="transparent",
            na.rm = F,
            show.legend=T) +
    scale_fill_manual(
                      values= natura_colors,
                       guide = guide_legend(
                                            override.aes = list(
                                                                linetype="solid",
                                                                shape = NA)
                                            ),
                       name="Natura2000"
                       )+
    theme_bw()

ggsave("../figures/map_natura.png", 
           plot=g_base_n2000, 
           height = 20, 
           width = 20,
           dpi = 300, 
           units="cm",
           device="png")

### natura2000 with all points
g_art17_n2000 <- g_base_n2000 +
    geom_sf(species_samples_presence_final_sf,
            mapping=aes(
                        shape=basisOfRecord,
                        color=datasetName),
            size=2,
            alpha=0.8,
            show.legend=T) +
    scale_color_manual(values=datasets_colors,
                        name = "Datasets")+
    guides(
           fill=guide_legend(position = "right",override.aes = list(linetype = 0,color=NA)),
           shape=guide_legend(position = "right"),
           color=guide_legend(position = "right",override.aes = list(linetype = 0,fill=NA)))+
    theme(legend.position = "right")


ggsave("../figures/map_art17_invertebrates_natura.png", 
           plot=g_art17_n2000, 
           height = 30, 
           width = 40,
           dpi = 300, 
           units="cm",
           device="png")

#----------------------------------------------------------------------------#
############################## Species Occurrences ##############################
#----------------------------------------------------------------------------#

### figures of each invertebrate of art17 for Greece

for (i in seq_along(species_with_data)){
    species_occurrences <- species_samples_presence_final_sf |>
        filter(species==species_with_data[i])

    dataset_colors_f <- datasets_colors[unique(species_occurrences$datasetName)]
    print(species_with_data[i])

    
    species_gr_map <- g_base_n2000 +
        geom_sf(species_occurrences,
                mapping=aes(
                        shape=basisOfRecord,
                        color=datasetName),
                size=2,
                alpha=0.9,
                show.legend=T) +
        scale_color_manual(values=dataset_colors_f,
                        name = "Datasets") +
        guides(
               fill=guide_legend(position = "right",override.aes = list(linetype = 0,color=NA)),
               shape=guide_legend(position = "right"),
               color=guide_legend(position = "right",override.aes = list(linetype = 0,fill=NA)))+
        ggtitle(paste(species_with_data[i]))+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position = "right",
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_occurrences.png", sep=""), 
           plot=species_gr_map, 
           height = 28, 
           width = 35,
           dpi = 300, 
           units="cm",
           device="png")

}

############################## Species Distribution ##############################

for (i in seq_along(species_with_data)){
    # use the 10X10 as a proxy of distribution

    # distribution 
    species_dist <- species_dist_national_rep_wgs_s |>
        filter(submittedName==species_with_data[i])

    # find the grids that don't have points
    dist_int <- st_intersects(species_dist, species_occurrences)
    no_points_index <- lengths(dist_int) == 0 
    points_index <- lengths(dist_int) > 0 
    polygons_without_points <- species_dist[no_points_index, ] 
    polygons__points <- species_dist[points_index, ] 

    
    # summary of points with grids
    locations_10_grid_samples <- st_join(eea_10km_wgs, species_occurrences, left=F) |>
        distinct(geometry,CELLCODE, decimalLatitude, decimalLongitude) |>
        group_by(geometry,CELLCODE) |>
        summarise(n_samples=n(),.groups="keep")

    # colors of datasets present
    dataset_colors_f <- datasets_colors[unique(species_occurrences$datasetName)]

    print(species_with_data[i])

    # Maps
    species_gr_map <- ggplot()+
        geom_sf(greece_regions, mapping=aes()) +
        geom_sf(N2000_v32_wgs, mapping=aes(fill=SITETYPE),
                alpha=0.6,
                #colour="transparent",
                na.rm = F,
                show.legend=T) +
        scale_fill_manual(
                          values= natura_colors,
                          guide = guide_legend(
                                                override.aes = list(alpha=1,
                                                                    linetype="solid",
                                                                    shape = NA)
                                                ),
                           name="Natura2000"
                           )+
        guides(
               fill=guide_legend(position = "inside",override.aes = list(alpha=1, linetype = 0,color=NA)),
               color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)))+
        new_scale_fill()+
        geom_sf(polygons_without_points, mapping=aes(fill="Cells without points"),
                alpha=0.8,
                colour="transparent",
                na.rm = F) +
        guides(
               fill=guide_legend(position = "inside",override.aes = list(linetype = 0,color=NA)))+
        new_scale_fill()+
        geom_sf(locations_10_grid_samples, mapping=aes(fill=n_samples),
                alpha=0.8,
                colour="transparent",
                na.rm = F) +
        scale_fill_gradient(low="gray50",
                            high="gray1",
                            guide = "colourbar")+
        guides(
               color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)),
               fill = guide_colourbar(position = "inside",
                                      alpha = 1,
                                      ticks = F,
                                      label = T,
                                      title="n_samples",
                                      title.vjust = 0.8,
                                      order = 1))+
        geom_point(species_occurrences,
                mapping=aes(x=decimalLongitude,
                            y=decimalLatitude,
                            color=datasetName),
                size=1,
                alpha=0.6) +
        scale_color_manual(values=dataset_colors_f,
                            name = "Datasets")+
        guides(
               fill=guide_legend(position = "inside",override.aes = list(linetype = 0,color=NA)))+
        ggtitle(paste(species_with_data[i]))+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position.inside = c(0.87, 0.65),
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_distribution.png", sep=""), 
           plot=species_gr_map, 
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")

}
############################## Species Population ##############################
### population is defined here as the number of 1X1km cells
### remove GBIF
########## Merge species occurrence data with previous distribution data#############

for (i in seq_along(species_with_data)){
    # use the 1X1 as a proxy of population

    species_occurrences <- species_samples_art17 |>
        filter(species==species_with_data[i]) |>
        filter(datasetName!="GBIF")
    
    dataset_colors_f <- datasets_colors[unique(species_occurrences$datasetName)]

    print(species_with_data[i])

    # Maps
    species_gr_map <- g_base_n2000+
        geom_sf(locations_1_grid_samples, mapping=aes(fill=n_samples),
                alpha=0.8,
                colour="transparent",
                na.rm = F) +
        scale_fill_gradient(low="gray50",
                            high="gray1",
                            guide = "colourbar")+
        guides(
               color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)),
               fill = guide_colourbar(position = "inside",
                                      alpha = 1,
                                      ticks = F,
                                      label = T,
                                      title="n_samples",
                                      title.vjust = 0.8,
                                      order = 1))+
        geom_point(species_occurrences,
                mapping=aes(x=decimalLongitude,
                            y=decimalLatitude,
                            color=datasetName),
                size=0.2,
                alpha=0.6) +
        scale_color_manual(values=dataset_colors_f,
                            name = "Datasets")+
        guides(
               fill=guide_legend(position = "inside",override.aes = list(linetype = 0,color=NA)))+
        ggtitle(paste(species_with_data[i]))+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position.inside = c(0.87, 0.65),
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_population.png", sep=""), 
           plot=species_gr_map, 
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")

}


################################ Species Range ##############################

# calculate range for all species

species_range = list()
polygons_no_points_all <- list()
colors_cell_origin <- c("range"="mediumvioletred",
                        "distribution"="limegreen",
                        "distribution orphan cell"="lightsalmon4")

for (i in seq_along(species_with_data)){
    #i=1
    species_occurrences <- species_samples_art17 |>
        filter(species==species_with_data[i]) |>
        filter(datasetName!="GBIF")
    
    # distribution 
    species_dist <- species_dist_national_rep_wgs_s |>
        filter(submittedName==species_with_data[i])

    # find the grids that don't have points
    dist_int <- st_intersects(species_dist, species_occurrences)
    no_points_index <- lengths(dist_int) == 0 
    points_index <- lengths(dist_int) > 0 
    polygons_without_points <- species_dist[no_points_index, ] 
    polygons__points <- species_dist[points_index, ] 
    polygons_no_points_all[[i]]<- polygons_without_points

    polygons_without_points_centroids <- st_centroid(polygons_without_points)

    eea_polygons_without_points <- st_intersects(eea_10km_wgs,polygons_without_points_centroids)
    eea_polygons_without_points_a <- eea_10km_wgs[lengths(eea_polygons_without_points) > 0,] |>
        mutate(cell_origin="distribution orphan cell")

    dataset_colors_f <- datasets_colors[unique(species_occurrences$datasetName)]
    locations_10_grid_samples <- st_join(eea_10km_wgs, species_occurrences, left=F) |>
        distinct(geometry,CELLCODE, decimalLatitude, decimalLongitude) |>
        group_by(geometry,CELLCODE) |>
        summarise(n_samples=n(),.groups="keep")
    
    grids <- eea_10km_wgs |>
        filter(CELLCODE %in% locations_10_grid_samples$CELLCODE) |>
        mutate(cell_origin="distribution") |>
        rbind(eea_polygons_without_points_a)
    
    expanded_range <- expand_range_with_gap_distance(distribution = grids,
                                                     full_grid = eea_10km_wgs,
                                                     gap_distance_m = 40000)

    species_range[[species_with_data[i]]] <- expanded_range


    cell_origin_colors_f <- colors_cell_origin[unique(expanded_range$cell_origin)]

    # Maps
    species_gr_map <- ggplot()+
        geom_sf(greece_regions, mapping=aes()) +
        geom_sf(N2000_v32_wgs, mapping=aes(fill=SITETYPE),
                alpha=0.6,
                #colour="transparent",
                na.rm = F,
                show.legend=T) +
        scale_fill_manual(
                          values= natura_colors,
                          guide = guide_legend(
                                                override.aes = list(alpha=1,
                                                                    linetype="solid",
                                                                    shape = NA)
                                                ),
                           name="Natura2000"
                           )+
        guides(
               fill=guide_legend(position = "inside",override.aes = list(alpha=1, linetype = 0,color=NA)),
               color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)))+
        new_scale_fill()+
        geom_sf(polygons_without_points, mapping=aes(fill="Cells without points"),
                alpha=0.8,
                colour="transparent",
                na.rm = F) +
        guides(
               fill=guide_legend(position = "inside",override.aes = list(linetype = 0,color=NA)))+
        new_scale_fill()+
        geom_sf(expanded_range, mapping=aes(fill=cell_origin),
                alpha=0.8,
                colour="transparent",
                na.rm = F) +
        scale_fill_manual(
                          values=cell_origin_colors_f ,
                          guide = guide_legend(
                                                override.aes = list(
                                                                    linetype="solid",
                                                                    shape = NA)
                                                ),
                           name="Cell origin"
                           )+
        geom_point(species_occurrences,
                mapping=aes(x=decimalLongitude,
                            y=decimalLatitude,
                            color=datasetName),
                size=1,
                alpha=0.6,
                show.legend=T) +
        scale_color_manual(values=dataset_colors_f,
                            name = "Datasets")+
        guides(
               fill=guide_legend(position = "inside",override.aes = list(linetype = 0,color=NA)))+
        ggtitle(paste(species_with_data[i]))+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position.inside = c(0.87, 0.75),
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_range.png", sep=""), 
           plot=species_gr_map, 
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")
    
}

