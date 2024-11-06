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
library(units)
library(vegan)


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
