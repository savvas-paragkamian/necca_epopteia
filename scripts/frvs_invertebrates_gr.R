#!/usr/bin/env Rscript

## Script name: frvs_invertebrates_gr.R
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

gr_1km <- sf::st_read("../spatial_data/eea_reference_grid/gr_1km.shp") |>
    st_transform(., crs="WGS84")

gr_10km <- sf::st_read("../spatial_data/eea_reference_grid/gr_10km.shp") |>
    st_transform(., crs="WGS84")


N2000_v32 <- sf::st_read("../spatial_data/N2000_spatial_GR_2021_12_09_v32/N2000_spatial_GR_2021_12_09_v32.shp")

N2000_v32_wgs <- st_transform(N2000_v32,4326)

########################### Load Species Data ###########################
### Species occurrences enriched ######
species_occurrences_invertebrates <- read_delim("../results/species_occurrences_invertebrates.tsv",delim="\t")
species_occurrences_spatial <- read_delim("../results/species_occurrences_spatial.tsv",delim="\t")
species_samples_art17 <- read_delim("../results/species_samples_art17.tsv", delim="\t")

species_taxonomy <- read_delim("../results/species_gbif_taxonomy_curated.tsv",delim="\t")
iucn_art17_invert_all <- read_delim("../results/iucn_art17_invert_all.tsv", delim="\t")
iucn_art17_invert_no_tax <- iucn_art17_invert_all |>
    dplyr::select(-ends_with("Name"),scientificName)

sspecies_dist_national_rep <- sf::st_read("../spatial_data/National report_2013_2018_shp/GR_Art17_species_distribution.shp")
species_dist_national_rep_sens <- sf::st_read("../spatial_data/National report_2013_2018_shp/GR_Art17_species_distribution_sensitive.shp")


########################## Flowchart for FRVs ##########################
# keep only one species name from the synonyms
species_info <- species_samples_art17 |>
    distinct(species) |> 
    rename("verbatim_name"="species") |>
    left_join(species_taxonomy, by=c("verbatim_name"="verbatim_name")) |>
    left_join(iucn_art17_invert_no_tax, by=c("species"="scientificName"))
# species spatial 
species_art17_spatial <- species_occurrences_spatial |>
    rename("verbatim_name"="species") |>
    filter(verbatim_name %in% species_taxonomy$verbatim_name) |>
    left_join(species_info, by=c("verbatim_name"="verbatim_name")) |>
    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),
             remove=F,
             crs="WGS84")


species_art17_spatial_f <- species_art17_spatial |>
    dplyr::select(verbatim_name,basisOfRecord,datasetName)

dir.create("../results/species_art17_spatial", recursive = TRUE, showWarnings = FALSE)
st_write(species_art17_spatial_f,"../results/species_art17_spatial/species_art17_spatial.shp",delete_layer = TRUE)

# summary of resourses
#
resources_summary_art17 <- species_samples_art17 |>
    group_by(datasetName,species) |>
    summarise(n_occurrences=n(), .groups="keep") |>
    group_by(datasetName) |>
    summarise(n_occurrences=sum(n_occurrences), n_species=n())

# what is known for each species

species_points <- species_art17_spatial |>
    distinct(species,decimalLatitude,decimalLongitude) |>
    group_by(species) |>
    summarise(n_points=n())

species_1km <- species_art17_spatial |>
    distinct(species,CELLCODE_eea_1km) |>
    group_by(species) |>
    summarise(n_1km=n())

species_1km_n2000 <- species_art17_spatial |>
    distinct(species,CELLCODE_eea_1km,SITECODE_N2000_v32_scispa,SITECODE_N2000_v32_spa,SITECODE_N2000_v32_sci) |>
    group_by(species) |>
    summarise(across(starts_with("SITECODE"), ~sum(!is.na(.)), .names = "count_{.col}"))


species_summary <- species_info |>
    dplyr::select(-verbatim_name) |>
    distinct() |>
    left_join(species_points) |>
    left_join(species_1km) |>
    left_join(species_1km_n2000)


write_delim(species_summary, "../results/species_summary.tsv", delim="\t")

########################## Parnassius apollo ###########################


parnassius_dist <- sf::st_read("../data/Parnassius apollo AP 2019/AP_Papollo_Distribution_LAEA.shp") |>
    st_transform(crs="WGS84")

species_samples_art17_parnasious <- species_art17_spatial |>
    filter(species=="Parnassius apollo") |>
    filter(X_eudem_dem_4258_europe>600)

apollo_mean <- species_samples_art17_parnasious |>
    dplyr::select(-starts_with("X_hilda_plu")) |>
    summarise(
    across(
      where(is.numeric),
      list(mean = mean),
      na.rm = TRUE
    )
  )


write_delim(apollo_mean,"../results/apollo_mean.tsv", delim="\t")

### hotspot
locations_10_grid_samples <- st_join(gr_10km, species_samples_art17_parnasious, left=F) |>
    distinct(geometry,CELLCODE, decimalLatitude, decimalLongitude) |>
    group_by(geometry,CELLCODE) |>
    summarise(n_samples=n(),.groups="keep")


hotspot_apollo_map <- ggplot() +
    geom_sf(greece_regions, mapping=aes()) +
    geom_sf(N2000_v32_wgs, mapping=aes()) +
    geom_sf(locations_10_grid_samples, mapping=aes(fill=n_samples),
            alpha=0.8,
            colour="transparent",
            na.rm = F,
            show.legend=T) +
    
    geom_sf(species_samples_art17_parnasious, mapping=aes(),size=1,alpha=0.4) +
    scale_fill_gradient(low="#f0e442",
                        high="#d55e00",
                        guide = "colourbar")+
    coord_sf(crs="wgs84") +
    guides(fill = guide_colourbar(ticks = F,
                                  label = T,
                                  title="n_samples",
                                  title.vjust = 0.8,
                                  order = 1))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../figures/hotspots_parnassius_apollo_map.png", 
       plot=hotspot_apollo_map, 
       height = 30, 
       width = 30,
       dpi = 200, 
       units="cm",
       device="png")
