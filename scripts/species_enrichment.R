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

invertebrates_92_43 <- readxl::read_excel("../data/invertebrates_92_43.xlsx", skip=2, col_names=F)

colnames(invertebrates_92_43) <- c("no","area","SPECIES_NAME","SPECIES_ID","ORDER","CLASS","PRIORITY","ANNEX_II","ANNEX_IV","ANNEX_V","KD","POPULATION_TREND","POPULATION_SIZE_UNIT","OCCURRENCE","SD")

species_names <- unique(invertebrates_92_43$SPECIES_NAME) 

#species_names <- c("Apatura metis",
#                   "Astacus astacus",
#                   "Austropotamobius torrentium",
#                   "Bolbelasmus unicornis",
#                   "Buprestis splendens",
#                   "Catopta thrips",
#                   "Cerambyx cerdo",
#                   "Coenagrion ornatum",
#                   "Cordulegaster heros",
#                   "Dioszeghyana schmidtii",
#                   "Eriogaster catax",
#                   "Euphydryas aurinia",
#                   "Euplagia quadripunctaria",
#                   "Hyles hippophaes",
#                   "Lindenia tetraphylla",
#                   "Lucanus cervus",
#                   "Lycaena dispar",
#                   "Maculinea arion",
#                   "Ophiogomphus cecilia",
#                   "Papilio alexanor",
#                   "Paracaloptenus caloptenoides",
#                   "Parnassius apollo",
#                   "Parnassius mnemosyne",
#                   "Polyommatus eroides",
#                   "Probaticus subrugosus",
#                   "Proserpinus proserpina",
#                   "Pseudophilotes bavius",
#                   "Rhysodes sulcatus",
#                   "Rosalia alpina",
#                   "Stenobothrus eurasius",
#                   "Stylurus flavipes",
#                   "Unio crassus",
#                   "Unio elongatulus",
#                   "Vertigo angustior",
#                   "Vertigo moulinsiana",
#                   "Zerynthia polyxena",
#                   "Morimus asper funereus",
#                   "Osmoderma eremita Complex",
#                   "Hirudo verbana")

## Natura2000 v32 version

natura_v32_species <- read_delim("../data/Natura2000DB_V32_species.tsv",delim="\t")
natura_v32_site <- read_delim("../data/Natura2000DB_V32_site.tsv",delim="\t")
natura_v32_region <- read_delim("../data/Natura2000DB_V32_region.tsv",delim="\t")
natura_v32_other_species <- read_delim("../data/Natura2000DB_V32_other_species.tsv",delim="\t")

groups_df <- data.frame(SPECIES_GROUP=c("R","B","A","M","I","F","P"),
                    groups_names=c("Reptiles", "Birds", "Amphibians", "Mammals", "Invertebrates", "Fish", "Plants"))
site_types_df <- data.frame(SITE_TYPE=c("A","B","C"), SITE_TYPE_NAME=c("SPAs","SCIs_SACs", "SPAs_and_SCIs_SACs"))

################################ Descriptives ###################################
##  
all_species <- unique(c(unique(natura_v32_species$SPECIES_NAME),unique(natura_v32_other_species$OTHER_SPECIES_NAME)))

natura_v32_species_sum <- natura_v32_species |> 
    distinct(SPECIES_NAME,SPECIES_GROUP,SITE_CODE) |>

natura_v32_other_species_sum <- natura_v32_other_species |> 
    distinct(OTHER_SPECIES_NAME,OTHER_SPECIES_GROUP, SITE_CODE)

colnames(natura_v32_other_species_sum) <- colnames(natura_v32_species_sum)

natura_v32_all_species <- rbind(natura_v32_other_species_sum, natura_v32_species_sum) |>
    distinct()

### species
groups_species_summary <- natura_v32_all_species |>
    distinct(SPECIES_NAME, SPECIES_GROUP) |>
    group_by(SPECIES_GROUP) |> 
    summarise(n_species=n()) |>
    left_join(groups_df)

invertebrates_all <- natura_v32_all_species |>
    filter(SPECIES_GROUP=="I")

invertebrates_all_natura_summary <- invertebrates_all |>
    group_by(SITE_CODE) |>
    summarise(n_species=n()) |> 
    left_join(natura_v32_region)

invertebrates_all_natura_summary_region <- invertebrates_all_natura_summary |>
    group_by(REGION_NAME) |>
    summarise(n_species=sum(n_species), n_sites=n())

### regions
## some sites are in multiple regions 
# > natura_v32_region |> distinct(REGION_CODE,SITE_CODE) |> group_by(SITE_CODE) |> summarise(n=n()) |> arrange(desc(n))

## overlap
## 16 not included, 23 included
invertebrates_not_in_natura <- species_names[which(!(species_names %in% natura_v32_species$SPECIES_NAME))]
invertebrates_not_in_other_natura <- species_names[which(!(species_names %in% natura_v32_other_species$OTHER_SPECIES_NAME))]


invertebrates_not_everywhere <- species_names[which(!(species_names %in% all_species))] 

species_names %in% unique(natura_v32_other_species$OTHER_SPECIES_NAME)

## Species from excel statistics

natura_v32_all_species_excel <- natura_v32_all_species |>
    filter(SPECIES_NAME %in% species_names) |>
    group_by(SPECIES_NAME) |>
    summarise(n_sites=n(),SITE_CODE=str_c(SITE_CODE, collapse = ","))

write_delim(natura_v32_all_species_excel, "../results/natura_v32_all_species_excel.tsv", delim="\t")

### sites

groups_sites_summary <- natura_v32_all_species |>
    distinct(SITE_CODE, SPECIES_GROUP,SPECIES_NAME) |>
    group_by(SITE_CODE,SPECIES_GROUP) |> 
    summarise(n_species=n(),
              species_names=str_c(SPECIES_NAME,collapse = ","),
              .groups="keep") |>
    left_join(natura_v32_site) |>
    left_join(site_types_df)

write_delim(groups_sites_summary,
            "../results/natura_groups_summary_all.tsv",
            delim="\t")

groups_sites_summary_invertebrates <- groups_sites_summary |>
    filter(SPECIES_GROUP=="I")

write_delim(groups_sites_summary_invertebrates,
            "../results/natura_groups_summary_invertebrates.tsv",
            delim="\t")

#### sites
groups_site_type_summary <- natura_v32_all_species |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME, SPECIES_GROUP,SITE_CODE) |>
    left_join(groups_df) |>
    group_by(SITE_TYPE_NAME,groups_names) |>
    summarise(n_sites=n(), .groups="keep") 

groups_site_type_sites_summary_all <- natura_v32_all_species |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME, SITE_CODE) |>
    group_by(SITE_TYPE_NAME) |>
    summarise(n_sites=n()) |>
    mutate(groups_names="All") |>
    rbind(groups_site_type_summary)

#### type species
groups_site_type_species_summary <- natura_v32_all_species |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME, SPECIES_GROUP,SPECIES_NAME) |>
    left_join(groups_df) |>
    group_by(SITE_TYPE_NAME,groups_names) |>
    summarise(n_species=n(), .groups="keep") 

groups_site_type_species_summary_all <- natura_v32_all_species |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME,SPECIES_NAME) |>
    group_by(SITE_TYPE_NAME) |>
    summarise(n_species=n()) |>
    mutate(groups_names="All") |>
    rbind(groups_site_type_species_summary)


### species from annex II
groups_sites_invertebrates <- groups_sites_summary |>
    filter(SPECIES_GROUP=="I")

species_annex_II <- invertebrates_92_43 |>
    filter(!is.na(ANNEX_II))

groups_site_type_species_annex_II <- natura_v32_all_species |>
    filter(SPECIES_NAME %in% species_annex_II$SPECIES_NAME) |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME, SPECIES_GROUP,SPECIES_NAME) |>
    mutate(groups_names="Invertebrates_annex_II") |>
    group_by(SITE_TYPE_NAME,groups_names) |>
    summarise(n_species=n(), .groups="keep") 

groups_site_type_annex_II_summary <- natura_v32_all_species |>
    filter(SPECIES_NAME %in% species_annex_II$SPECIES_NAME) |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME, SPECIES_GROUP,SITE_CODE) |>
    mutate(groups_names="Invertebrates_annex_II") |>
    group_by(SITE_TYPE_NAME,groups_names) |>
    summarise(n_sites=n(), .groups="keep") 

natura_sites_groups_species_annex_II <- groups_site_type_species_annex_II |>
    left_join(groups_site_type_annex_II_summary,
    by=c("SITE_TYPE_NAME","groups_names"))

### sites and species together

natura_sites_groups_species_summary <- groups_site_type_sites_summary_all |>
    left_join(groups_site_type_species_summary_all,
              by=c("SITE_TYPE_NAME","groups_names")) |>
    rbind(natura_sites_groups_species_annex_II)

write_delim(natura_sites_groups_species_summary,
            "../results/natura_sites_groups_species_summary.tsv",
            delim="\t")

################################ Enrichment ###################################
# GBIF retrieve data for all arthropod species that have been assessed in IUCN
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
