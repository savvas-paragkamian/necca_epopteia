#!/usr/bin/Rscript

## Script name: species_enrichment.R
##
## Purpose of script: use public databases to retrieve information of 
## species regarding their taxonomy, their global distribution and their 
## IUCN status
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-11-06

library(sf)
library(tidyverse)
library(readxl)
library(taxize)
library(rgbif)
library(units)
library(vegan)
library(rnaturalearth)

# Load data
## borders for maps
world <- ne_countries(scale = "medium", returnclass = "sf")

country_name <- "Greece"
greece <- world[world$name == country_name, ]

hellenic_borders_shp <- sf::st_read("../spatial_data/hellenic_borders/hellenic_borders.shp") |>
    st_transform(4326)

# the species of the previous epopteia 2015
species_names <- c("Apatura metis",
                   "Astacus astacus",
                   "Austropotamobius torrentium",
                   "Bolbelasmus unicornis",
                   "Buprestis splendens",
                   "Catopta thrips",
                   "Cerambyx cerdo",
                   "Coenagrion ornatum",
                   "Cordulegaster heros",
                   "Dioszeghyana schmidtii",
                   "Eriogaster catax",
                   "Euphydryas aurinia",
                   "Euplagia quadripunctaria",
                   "Hyles hippophaes",
                   "Lindenia tetraphylla",
                   "Lucanus cervus",
                   "Lycaena dispar",
                   "Maculinea arion",
                   "Ophiogomphus cecilia",
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
                   "Stenobothrus eurasius",
                   "Stylurus flavipes",
                   "Unio crassus",
                   "Unio elongatulus",
                   "Vertigo angustior",
                   "Vertigo moulinsiana",
                   "Zerynthia polyxena",
                   "Morimus asper funereus",
                   "Osmoderma eremita Complex",
                   "Hirudo verbana")

## GBIF retrieve data for all arthropod species that have been assessed in IUCN
### NOT run takes time. 
gbif_species <- get_gbifid(species_names,ask=F)


species_gbif_df <- data.frame(sci_name=species_names, gbifid=gbif_species)
write_delim(species_gbif_df, "../results/species_gbif.tsv", delim="\t")
#### takes even more time!!!
classification_species <- classification(species_gbif_df$gbifid.ids, db = 'gbif')

classification_species_d <- do.call(rbind, classification_species) |>
    rownames_to_column(var="gbif") |> 
    mutate(gbif = gsub("\\.(.*)","", gbif)) |>
    dplyr::select(-id) |>
    distinct() |>
    na.omit(gbif) |>
    pivot_wider(id_cols=gbif, names_from=rank, values_from=name ) |>
    mutate(gbif=as.numeric(gbif)) 

write_delim(classification_species_d, "../results/classification_species_gbif.tsv", delim="\t")

# Resolve names
## gnr_datasources() |> filter(title=="GBIF Backbone Taxonomy") id=11
gnr_species <- gnr_resolve(species_gbif_df$sci_name)
gnr_species_gbif <- gnr_resolve(species_gbif_df$sci_name, data_source_ids=11)

write_delim(gnr_species, "../results/gnr_species_names.tsv", delim="\t")

# GBIF occurrences
## need to set GBIF credential in the .Renviron file

gbif_taxon_keys <- as.numeric(na.omit(species_gbif_df$gbifid.ids))

## run once to request the download from the server
#occ_download(
#pred_in("taxonKey", gbif_taxon_keys), # important to use pred_in
#pred("hasCoordinate", TRUE),
#pred("hasGeospatialIssue", FALSE),
#format = "SIMPLE_CSV"
#)

# to check the status of the download
# occ_download_wait('0026745-241024112534372') 

gbif_species_occ <- occ_download_get('0026745-241024112534372') |>
    occ_download_import()

write_delim(gbif_species_occ, "../results/gbif_species_occ.tsv", delim="\t")


gbif_species_occ_sf <- gbif_species_occ |> 
    st_as_sf(coords=c("decimalLatitude","decimalLongitude"),
             remove=F,
             crs="WGS84") |>
    st_transform(4326)

hellenic_borders_bbox <- st_as_sf(st_as_sfc(st_bbox(hellenic_borders_shp)))

gbif_species_world_map <- ggplot() +
    geom_sf(world, mapping=aes()) +
    geom_point(gbif_species_occ_sf,
            mapping=aes(x=decimalLongitude,
                        y=decimalLatitude,
                        color=species),
            size=1.8,
            alpha=0.8,
            show.legend=T) +
    geom_sf(hellenic_borders_shp, mapping=aes()) +
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../figures/gbif_species_world_map.png", 
       plot=gbif_species_world_map, 
       height = 40, 
       width = 60,
       dpi = 300, 
       units="cm",
       device="png")

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
    geom_sf(hellenic_borders_shp, mapping=aes()) +
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
