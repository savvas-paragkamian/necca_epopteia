#!/usr/bin/env Rscript

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
library(rnaturalearth)

############################# Load data ############################
## borders for maps
greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")

species_taxonomy <- read_delim("../results/species_gbif_taxonomy_curated.tsv",delim="\t")

species_names_combined <- as.character(species_taxonomy$verbatim_name)

################################ Enrichment ###################################

####################### edaphobase ########################
edaphobase_gr <- read_delim("../data/2025-02-26-edaphobase-export_GR.csv", delim=";")

# searching in the df multiple tokens using grepl
edaphobase_gr_art17 <- edaphobase_gr[Reduce(`|`, lapply(species_names_combined, function(p) grepl(p, edaphobase_gr$`Valid taxon`,ignore.case = T))), ]


############################################################################################
#################### Species Occurrences Data Homogenisation ###############################
############################################################################################
##### Gbif data
## data filter for coordinate precision
print("gbif")
##
gbif_all <- read_delim("../data/gbif_invertebrate_species_occ.tsv", delim="\t") 

gbif_species_occ <- read_delim("../data/gbif_invertebrate_species_occ.tsv", delim="\t") |>
    mutate(datasetName="GBIF") |>
    mutate(
           species = ifelse(!is.na(verbatimScientificName) &
                                    verbatimScientificName=="Panaxia quadripunctaria",
                                     "Euplagia quadripunctaria",
                                     as.character(species))
    ) |>
    rename("submittedName"="species") |>
    filter(coordinateUncertaintyInMeters<1000)


## Define the bounding box coordinates
#xmin <- 19.37359
#ymin <- 34.80202
#xmax <- 29.64306
#ymax <- 41.7485

gbif_species_occ_sf <- gbif_species_occ |>
    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),
             remove=F,
             crs="WGS84")


within_mat <- st_intersects(gbif_species_occ_sf,greece_regions, sparse = FALSE)

gbif_species_occ_sf <- gbif_species_occ_sf[rowSums(within_mat) > 0, ]

gbif_species_occ_gr <- gbif_species_occ_sf |>
    st_drop_geometry()

print("end gbif")

######## NECCA compilation of previous Monitoring Data not included in ENVECO database
E1X_MDPP_2014_2024_samples_data <- read_xlsx("../data/Ε1Χ_ΒΔ_ΠΡΩΤΟΓΕΝΩΝ_ΦΔ+ΜΔΠΠ_2014-2024.xlsx",
                           sheet="Δείγματα Ασπόνδυλων",
                           col_names=T) |> slice(-1) |> 
    mutate(decimalLatitude=as.numeric(`Γεωγραφικό Πλάτος (WGS84) Αρχη`),
           decimalLongitude=as.numeric(`Γεωγραφικό Μήκος (WGS84) Αρχή`)) |> 
    filter(!is.na(decimalLongitude)) |>
    mutate(datasetName="E1X_MDPP_2014_2024") |>
    mutate(basisOfRecord="MATERIAL_SAMPLE")


E1X_MDPP_2014_2024_species_data <- read_xlsx("../data/Ε1Χ_ΒΔ_ΠΡΩΤΟΓΕΝΩΝ_ΦΔ+ΜΔΠΠ_2014-2024.xlsx",
                           sheet="Είδη",
                           col_names=T
                           ) |> slice(-1)

E1X_MDPP_2014_2024_all <- E1X_MDPP_2014_2024_species_data |>
    mutate(species=if_else(`Όνομα είδους`=="Άλλο",
                           `Άλλο είδος`,
                           `Όνομα είδους`)) |>
    mutate(art17_92_43_EEC=if_else(`Όνομα είδους`!="Άλλο",
                                   TRUE,
                                   FALSE)) |>
    mutate(individualCount=as.numeric(`Αριθμός ατόμων είδους`)) |>
    left_join(E1X_MDPP_2014_2024_samples_data, by=c("Sam_ID"="Sam_ID")) |>
    mutate(submittedName=`Όνομα είδους`)

######### previous monitoring from ENVECO
##### references
E1X_DB_ref_samples_data <- read_xlsx("../data/Ε1Χ_ΒΔ_ΒΙΒΛΙΟΓΡΑΦΙΑΣ_ΑΣΠ_20241204.xlsx",
                                    sheet="Εξάπλωση ειδών και τ.ο.",
                                    col_names=T) |> slice(-1) |> filter(!is.na(`Κωδικός Αναφοράς`))

E1X_DB_refs_data <- read_xlsx("../data/Ε1Χ_ΒΔ_ΒΙΒΛΙΟΓΡΑΦΙΑΣ_ΑΣΠ_20241204.xlsx",
                                    sheet="Βιβλιογραφία",
                                    col_names=T) |> slice(-1) |> filter(!is.na(`Κωδικός Αναφοράς`))


E1X_DB_ref_samples_data$decimalLatitude <- as.numeric(E1X_DB_ref_samples_data$`Γεωγραφικό Πλάτος (WGS84)`)
E1X_DB_ref_samples_data$decimalLongitude <- as.numeric(E1X_DB_ref_samples_data$`Γεωγραφικό Μήκος (WGS84)`)
E1X_DB_ref_samples_data$species <- E1X_DB_ref_samples_data$`Ονομασία είδους`

E1X_DB_ref_all <- E1X_DB_ref_samples_data |>
    left_join(E1X_DB_refs_data,
              by=c("Κωδικός Αναφοράς"="Κωδικός Αναφοράς")) |>
    mutate(datasetName="E1X_DB_references") |>
    mutate(basisOfRecord="MaterialCitation") |>
    mutate(submittedName=`Ονομασία είδους`) |>
    mutate(individualCount=as.numeric(`Πλήθος ατόμων`))


##### samplings
E1X_DB_samples_data <- read_xlsx("../data/Ε1Χ_ΒΔ_ΠΡΩΤΟΓΕΝΩΝ_ΥΠ4_ΑΣΠΟΝΔΥΛΑ_20241204.xlsx",
                                    sheet="Δείγματα Ασπόνδυλων",
                                    col_names=T) |> slice(-1)

E1X_DB_samples_data$decimalLatitude <- as.numeric(E1X_DB_samples_data$`Γεωγραφικό Πλάτος (WGS84) Αρχη`)
E1X_DB_samples_data$decimalLongitude <- as.numeric(E1X_DB_samples_data$`Γεωγραφικό Μήκος (WGS84) Αρχή`)

E1X_DB_species_data <- read_xlsx("../data/Ε1Χ_ΒΔ_ΠΡΩΤΟΓΕΝΩΝ_ΥΠ4_ΑΣΠΟΝΔΥΛΑ_20241204.xlsx",
                                    sheet="Είδη",
                                    col_names=T) |> slice(-1) 

E1X_DB_all <- E1X_DB_species_data |> 
    filter(!is.na(Sam_ID)) |> 
    left_join(E1X_DB_samples_data, by=c("Sam_ID"="Sam_ID")) |>
    mutate(individualCount=as.numeric(`Αριθμός ατόμων είδους`)) |>
    mutate(submittedName=`Όνομα είδους`) |>
    mutate(datasetName = "E1X_DB") |>
    mutate(basisOfRecord="MATERIAL_SAMPLE")


## these data have all species that the researchers were looking for
## but they have found the species that have individualCount more than 0!
## so we are filtering them.
##
E1X_DB_select <- E1X_DB_all |>
    filter(individualCount>0)

######################## other private data

Invertebrates_records_Olga <- read_delim("../data/Invertebrates_records_Olga_20250427.csv", delim=",") |>
    mutate(decimalLongitude=Longitude,
           decimalLatitude=Latitude,
           submittedName=Species,
           individualCount=as.numeric(Individuals)) |>
    mutate(datasetName = "Invertebrates_records_Olga") |>
    mutate(basisOfRecord="MATERIAL_SAMPLE")

########################### NECCA Redlist ###########################

necca_redlist_points <- st_read("../data/necca_redlist/points_invertebrates.gpkg")

necca_redlist_points_df <- necca_redlist_points |>
    st_cast("POINT") %>% 
    mutate(decimalLongitude = st_coordinates(.)[, 1],
           decimalLatitude = st_coordinates(.)[, 2]) |>
    st_drop_geometry() |>
    rename("submittedName"="sci_name") |>
    mutate(datasetName = "NECCA_redlist") |>
    mutate(individualCount=NA) |>
    mutate(basisOfRecord="MATERIAL_SAMPLE")

necca_redlist_polygons <- st_read("../data/necca_redlist/polygons_invertebrates.gpkg")


########################### IUCN Redlist ###########################
#
iucn_summary <- read_delim("../data/redlist_species_data_c7359df6-ef0b-468c-b295-32816a14f56b/simple_summary.csv", delim=",")
iucn_assessments <- read_delim("../data/redlist_species_data_c7359df6-ef0b-468c-b295-32816a14f56b/assessments.csv", delim=",")

iucn_assessments_latest <- iucn_assessments |>
    group_by(scientificName) |>
    filter(yearPublished == max(yearPublished)) |>
    ungroup()

iucn_threats <- read_delim("../data/redlist_species_data_c7359df6-ef0b-468c-b295-32816a14f56b/threats.csv", delim=",")

iucn_points <- read_delim("../data/redlist_species_data_dbd309b7-fb96-4f22-91c2-f184787ada27/points_data.csv", delim=",")
#
iucn_art17_invert_points <- iucn_points |>
    filter(str_detect(sci_name,str_c(gsub(" .*","",species_names_combined), collapse = "|")))

# there are 25 species of Art 17 for greek invertebrates in iucn
iucn_art17_invert_summary <- iucn_summary |>
    filter(str_detect(scientificName,
                      str_c(species_names_combined, collapse = "|"))) |>
    distinct(scientificName,
             phylumName,
             orderName,
             className,
             familyName,
             genusName
             ) 


iucn_assessments_art17 <- iucn_assessments |>
    filter(str_detect(scientificName,str_c(species_names_combined, collapse = "|"))) |>
    group_by(scientificName) |>
    summarise(
              assessmentIds=paste(assessmentId, collapse="|"),
              scopes=paste(scopes, collapse="|"),
              redlistCategory=paste(redlistCategory, collapse="|"),
              population=paste(population, collapse="|"),
              populationTrend=paste(populationTrend, collapse="|")
              ) |>
    ungroup()

iucn_art17_invert_threats <- iucn_threats |>
    filter(str_detect(scientificName,str_c(species_names_combined, collapse = "|"))) |>
    group_by(scientificName) |>
    summarise(
              stressCode=paste(stressCode, collapse="|"),
              stressName=paste(stressName, collapse="|")
              ) |>
    ungroup()

## combine
iucn_art17_invert_all <- iucn_art17_invert_summary |>
    left_join(iucn_assessments_art17,
              by=c("scientificName"="scientificName")
              ) |>
    left_join(iucn_art17_invert_threats,
              by=c("scientificName"="scientificName")
              ) |>
    mutate(datasetName = "IUCN_redlist") |>
    rename("submittedName"="scientificName")

write_delim(iucn_art17_invert_all, "../results/iucn_art17_invert_all.tsv", delim="\t")
########################### Stenobothrus eurasius ###########################

stenobothrus_eurasius <- read_delim("../data/stenobothrus_eurasius.tsv",delim="\t")

############################# Species data integration ###################
### common column names
columns_to_keep <- c("submittedName",
                     "decimalLatitude",
                     "decimalLongitude",
                     "datasetName",
                     "basisOfRecord",
                     "individualCount")

species_occurrences_invertebrates <- list(gbif_species_occ_gr,
                                    E1X_MDPP_2014_2024_all,
                                    E1X_DB_select,
                                    E1X_DB_ref_all,
                                    necca_redlist_points_df,
                                    stenobothrus_eurasius,
                                    Invertebrates_records_Olga) |>
    map(~ dplyr::select(.x, all_of(columns_to_keep))) |>
    bind_rows() 


write_delim(species_occurrences_invertebrates, "../results/species_occurrences_invertebrates.tsv",delim="\t")

## only species art17

species_samples_art17_all <- species_occurrences_invertebrates |>
    filter(submittedName %in% species_names_combined) |>
    left_join(species_taxonomy, by=c("submittedName"="verbatim_name"))

write_delim(species_samples_art17_all,"../results/species_samples_art17_all.tsv", delim="\t")

