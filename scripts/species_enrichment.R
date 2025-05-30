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
                                    Invertebrates_records_Olga) |>
    map(~ dplyr::select(.x, all_of(columns_to_keep))) |>
    bind_rows() 


write_delim(species_occurrences_invertebrates, "../results/species_occurrences_invertebrates.tsv",delim="\t")

species_occurrences_invertebrates_sf <- species_occurrences_invertebrates |>
    filter(!is.na(decimalLongitude)) |>
    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),
             remove=F,
             crs="WGS84")


species_samples_art17_all <- species_occurrences_invertebrates |>
    filter(submittedName %in% species_names_combined) |>
    left_join(species_taxonomy, by=c("submittedName"="verbatim_name"))

species_samples_art17_sf <- species_samples_art17_all |>
    filter(!is.na(decimalLongitude)) |>
    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),
             remove=F,
             crs="WGS84")

### find the points that are on land and inside borders of greece
points <- species_samples_art17_sf

polygons <- greece_regions
intersects_mat <- st_intersects(points, polygons, sparse = FALSE)
# Use st_within to get logical vector of which points are inside the multipolygon
points_outside <- points[rowSums(intersects_mat) == 0, ]
points_inside_or_touching <- points[rowSums(intersects_mat) > 0, ]

# 1000 meters = 1 km
nearby_mat <- st_is_within_distance(points, polygons, dist = 2000, sparse = FALSE)

# Keep only points within 1 km of any polygon
points_nearby <- points[rowSums(nearby_mat) > 0, ]
points_away_2km <- points[rowSums(nearby_mat) == 0, ]



# plot

plot(st_geometry(polygons), col = NA, border = "blue")
plot(st_geometry(points), add = TRUE, col = "black", pch = 1)
plot(st_geometry(points_outside), add = TRUE, col = "red", pch = 16)
plot(st_geometry(points_nearby), add = TRUE, col = "yellow", pch = 16)
#plot(st_geometry(species_samples_art17_sf), add = TRUE, col = "blue", pch = 16)

species_samples_art17

write_delim(species_samples_art17,"../results/species_samples_art17.tsv", delim="\t")

#### save
species_samples_art17_open <- species_samples_art17 |>
    filter(datasetName!="Invertebrates_records_Olga")

write_delim(species_samples_art17_open,"../results/species_samples_art17_open.tsv", delim="\t")

species_with_data <- unique(species_samples_art17$species)

datasets_colors <- c(
                     "GBIF"="seagreen",
                     "NECCA_redlist"="#B31319",
                     "E1X_MDPP_2014_2024"="#FDF79C",
                     "E1X_DB"="#2BA09F",
                     "E1X_DB_references"="#141D43",
                     "Invertebrates_records_Olga"="#F85C29"
                     )
# load natura
N2000_v32 <- sf::st_read("../spatial_data/N2000_spatial_GR_2021_12_09_v32/N2000_spatial_GR_2021_12_09_v32.shp")

N2000_v32_wgs <- st_transform(N2000_v32,4326)


# base plot with Natura2000 areas of Greece
natura_colors <- c(
                   "SCI"="#E69F00",
                   "SPA"="#56B4E9",
                   "SCISPA"="#CC79A7"
)

g_base_n2000 <- ggplot()+
    geom_sf(greece_regions, mapping=aes()) +
    geom_sf(N2000_v32_wgs, mapping=aes(fill=SITETYPE),
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
    geom_point(species_samples_art17_open,
            mapping=aes(x=decimalLongitude,
                        y=decimalLatitude,
                        color=datasetName),
            size=1.2,
            alpha=0.8,
            show.legend=T) +
    scale_color_manual(values=datasets_colors,
                        name = "Datasets")+
    guides(
           fill=guide_legend(position = "inside",override.aes = list(linetype = 0,color=NA)),
           color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)))+
    theme(legend.position.inside = c(0.87, 0.75)
    )


ggsave("../figures/map_art17_invertebrates_natura.png", 
           plot=g_art17_n2000, 
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")

### figures of each invertebrate of art17 for Greece

for (i in seq_along(species_with_data)){
    species_occurrences <- species_samples_art17_sf |>
        filter(species==species_with_data[i])

    dataset_colors_f <- datasets_colors[unique(species_occurrences$datasetName)]
    print(species_with_data[i])

    
    species_gr_map <- g_base_n2000 +
        geom_point(species_occurrences,
                mapping=aes(x=decimalLongitude,
                            y=decimalLatitude,
                            color=datasetName),
                size=1.8,
                alpha=0.9,
                show.legend=T) +
        coord_sf(crs="WGS84") +
        scale_color_manual(values=dataset_colors_f,
                        name = "Datasets") +
        guides(
               fill=guide_legend(position = "inside",override.aes = list(linetype = 0,color=NA)),
               color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)))+
        ggtitle(paste(species_with_data[i]))+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position.inside = c(0.87, 0.75),
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_occurrences.png", sep=""), 
           plot=species_gr_map, 
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")

}


