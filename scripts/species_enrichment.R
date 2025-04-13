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

############################# Load data ############################
## borders for maps
greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")
hellenic_borders_shp <- sf::st_read("../spatial_data/hellenic_borders/hellenic_borders.shp")
########################## species information ##############################
# the species of the previous epopteia 2015

invertebrates_92_43 <- readxl::read_excel("../data/invertebrates_92_43.xlsx", skip=2, col_names=F)

colnames(invertebrates_92_43) <- c("no","area","SPECIES_NAME","SPECIES_ID","ORDER","CLASS","PRIORITY","ANNEX_II","ANNEX_IV","ANNEX_V","KD","POPULATION_TREND","POPULATION_SIZE_UNIT","OCCURRENCE","SD")

species_names <- unique(invertebrates_92_43$SPECIES_NAME) 

#### this list contains species names from multiple versions of Monitoring, EPOPTEIA I, EPOPTEIA II. 
species_names_combined <- c(
  "Apatura metis",
  "Astacus astacus",
  "Austropotamobius torrentium",
  "Bolbelasmus unicornis",
  "Buprestis splendens",
  "Callimorpha (Euplagia, Panaxia) quadripunctaria",
  "Catopta thrips",
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
  "Morimus funereus",
  "Ophiogomphus cecilia",
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
  "Vertigo angustior",
  "Vertigo moulinsiana",
  "Zerynthia polyxena"
)

################################ Enrichment ###################################

####################### edaphobase ########################
edaphobase_gr <- read_delim("../data/2025-02-26-edaphobase-export_GR.csv", delim=";")

# searching in the df multiple tokens using grepl
edaphobase_gr_art17 <- edaphobase_gr[Reduce(`|`, lapply(species_names_combined, function(p) grepl(p, edaphobase_gr$`Valid taxon`,ignore.case = T))), ]

######################## GBIF ########################
# GBIF retrieve data for all arthropod species that have been assessed in IUCN
### NOT run takes time. 
species_names <- as.character(invertebrates_all_species$SPECIES_NAME)
gbif_species <- get_gbifid(species_names,ask=F)


species_gbif_df <- data.frame(sci_name=species_names, gbifid=gbif_species)
write_delim(species_gbif_df, "../results/species_gbif_invertebrates.tsv", delim="\t")
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



############################################################################################
#################### Species Occurrences Data Homogenisation ###############################
############################################################################################
##### Gbif data

gbif_species_occ <- read_delim("../results/gbif_species_occ.tsv", delim="\t") |>
    mutate(datasetName="Gbif")

# Define the bounding box coordinates
xmin <- 19.37359
ymin <- 34.80202
xmax <- 29.64306
ymax <- 41.7485

gbif_species_occ_gr <- gbif_species_occ |>
    filter(decimalLongitude > 19.37359 & decimalLongitude<29.64306,
           decimalLatitude>34.80202 & decimalLatitude < 41.7485 )

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
    left_join(E1X_MDPP_2014_2024_samples_data, by=c("Sam_ID"="Sam_ID"))

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
    mutate(species=`Ονομασία είδους`) |>
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
    mutate(species=`Όνομα είδους`) |>
    mutate(datasetName = "E1X_DB") |>
    mutate(basisOfRecord="MATERIAL_SAMPLE")


## these data have all species that the researchers were looking for
## but they have found the species that have individualCount more than 0!
## so we are filtering them.
##
E1X_DB_select <- E1X_DB_all |>
    filter(individualCount>0)

######################## other private data

Invertebrates_records_Olga <- read_delim("../data/Invertebrates_records_Olga_20250316.csv", delim=",") |>
    mutate(decimalLongitude=Longitude,
           decimalLatitude=Latitude,
           species=Species,
           individualCount=as.numeric(Individuals)) |>
    mutate(datasetName = "Invertebrates_records_Olga") |>
    mutate(basisOfRecord="MATERIAL_SAMPLE")

########################### IUCN
#
#iucn <- read_delim("../data/redlist_species_data_dbd309b7-fb96-4f22-91c2-f184787ada27/points_data.csv", delim=",")
#
#iucn_92_43 <- iucn |>
#    filter(sci_name %in% species_names)
#
#iucn_parnassius <- iucn |>
#    filter(sci_name=="Parnassius apollo")


############################# Species data integration ###################
### common column names
columns_to_keep <- c("species",
                     "decimalLatitude",
                     "decimalLongitude",
                     "datasetName",
                     "basisOfRecord",
                     "individualCount")

species_samples_integration <- list(gbif_species_occ_gr,
                                    E1X_MDPP_2014_2024_all,
                                    E1X_DB_select,
                                    E1X_DB_ref_all,
                                    Invertebrates_records_Olga) |>
    map(~ dplyr::select(.x, all_of(columns_to_keep))) |>
    bind_rows()

write_delim(species_samples_integration, "../results/species_samples_integration.tsv",delim="\t")

species_samples_art17 <- species_samples_integration |>
    filter(species %in% species_names_combined)

species_samples_art17_sf <- species_samples_art17 |>
    filter(!is.na(decimalLongitude)) |>
    st_as_sf(coords=c("decimalLatitude","decimalLongitude"),
             remove=F,
             crs="WGS84")

species_with_data <- unique(species_samples_art17_sf$species)

for (i in seq_along(species_with_data)){
    species_occurrences <- species_samples_art17_sf |>
        filter(species==species_with_data[i])
    print(species_with_data[i])

    
    species_gr_map <- ggplot() +
        geom_sf(greece_regions, mapping=aes()) +
        geom_point(species_occurrences,
                mapping=aes(x=decimalLongitude,
                            y=decimalLatitude,
                            color=datasetName),
                size=1.8,
                alpha=0.8,
                show.legend=T) +
        coord_sf(crs="WGS84") +
        ggtitle(paste(species_with_data[i]))+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position = "bottom",
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/map_", species_with_data[i], "_occurrences.png", sep=""), 
           plot=species_gr_map, 
           height = 20, 
           width = 20,
           dpi = 300, 
           units="cm",
           device="png")

}






