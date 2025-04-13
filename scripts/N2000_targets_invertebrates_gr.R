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
library(units)
library(rnaturalearth)

############################# Load data ############################
## borders for maps
greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")

########################## species information ##############################
# the species of the previous epopteia 2015

invertebrates_92_43 <- readxl::read_excel("../data/invertebrates_92_43.xlsx", skip=2, col_names=F)

colnames(invertebrates_92_43) <- c("no","area","SPECIES_NAME","SPECIES_ID","ORDER","CLASS","PRIORITY","ANNEX_II","ANNEX_IV","ANNEX_V","KD","POPULATION_TREND","POPULATION_SIZE_UNIT","OCCURRENCE","SD")

species_names <- unique(invertebrates_92_43$SPECIES_NAME) 
##################### Natura2000 v32 version #############################

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
    distinct(SPECIES_NAME,SPECIES_GROUP,SITE_CODE) 

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

invertebrates_all_species <- invertebrates_all |>
    distinct(SPECIES_NAME)

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


