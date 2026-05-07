library(sf)
library(terra)
library(units)
library(readxl)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)


#### transform some sources 

build_national_report_distribution_grid_occurrences <- function(
  species_distribution,
  eea_grid_10km
) {

  species_dist_national_E1_eea <- sf::st_join(
    eea_grid_10km,
    species_distribution,
    join = sf::st_within,
    left = FALSE,
    largest = FALSE
  )

  centroids <- sf::st_centroid(species_dist_national_E1_eea)

  coords <- sf::st_coordinates(centroids)

  species_dist_national_E1_eea_coords <- cbind(
    species_dist_national_E1_eea,
    decimalLongitude = coords[, 1],
    decimalLatitude = coords[, 2]
  )

  species_dist_national_E1_eea <- species_dist_national_E1_eea_coords |>
    dplyr::distinct(
      submittedName,
      species,
      CELLCODE,
      decimalLongitude,
      decimalLatitude
    ) |>
    dplyr::rename(
      CELLCODE_eea_10km = CELLCODE
    ) |>
    dplyr::mutate(
      datasetName = "DistrMap_2013_2018",
      basisOfRecord = "ESTIMATED_CENTROID",
      recordNumber = paste0(
        "DistrMap_2013_2018_",
        sprintf("%04d", dplyr::row_number())
      ),
      collectionCode = "GR_Art17_species_distribution.shp",
      DistrMap_2013_2018 = TRUE
    )

  species_dist_national_E1_eea_coords_laea <- species_dist_national_E1_eea |>
    sf::st_as_sf(
      coords = c(
        "decimalLongitude",
        "decimalLatitude"
      ),
      remove = TRUE,
      crs = sf::st_crs(species_distribution)
    )

  species_dist_national_E1_eea_coords_wgs <- sf::st_transform(
    species_dist_national_E1_eea_coords_laea,
    4326
  )

  coords_wgs <- sf::st_coordinates(
    species_dist_national_E1_eea_coords_wgs
  )

  cbind(
    species_dist_national_E1_eea_coords_wgs,
    decimalLongitude = coords_wgs[, 1],
    decimalLatitude = coords_wgs[, 2]
  )
}

##### bind all data

################################## old code ###################

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
    mutate(individualCount = NA) 

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
eu_dem_gr <- rast("../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_3035_europe.tif")
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
    dplyr::select(species, datasetName,basisOfRecord,geometry,minimumDistanceFromBorders,eudem_dem_3035_europe) |>
    arrange(species)

write_delim(points_outside_cols,
            "../results/species_occurrences_outside_gr_land.tsv",
            delim="\t")


############################# Species data integration ###################
### common column names
columns_to_keep <- c("submittedName",
                     "decimalLatitude",
                     "decimalLongitude",
                     "datasetName",
                     "recordNumber",
                     "collectionCode",
                     "basisOfRecord",
                     "individualCount")

species_occurrences_invertebrates <- list(gbif_species_occ_gr,
                                    E1X_MDPP_2014_2024_all,
                                    E1X_DB_select,
                                    E1X_DB_ref_all,
                                    E2X_DB,
                                    necca_redlist_points_df,
                                    unio_crassus_complex_gr,
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


#### filtering






species_samples_simple <- species_samples_art17_all |> 
    filter(!is.na(decimalLatitude)) |>
    distinct(submittedName,
             decimalLatitude,
             decimalLongitude,
             datasetName,
             collectionCode,
             recordNumber,
             basisOfRecord,
             individualCount,
             species) |>
    group_by(across(-recordNumber)) |>
    summarise(recordNumber=str_c(recordNumber, collapse = ", "),.groups="keep") |>
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
    dplyr::select(-recordNumber) |>
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
             collectionCode,
             recordNumber,
             basisOfRecord,
             individualCount,
             CELLCODE_eea_10km,
             species)



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
                                       eudem_dem_3035_europe <= 450 & 
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


#species_samples_presence_final <- read_delim("../results/species_samples_presence_final.tsv", delim="\t")

################### with private data #################
species_samples_presence_final_private <- species_samples_presence_pop |>
    st_drop_geometry()

write_delim(species_samples_presence_final_private,
            "../results/species_samples_presence_final_private.tsv",
            delim="\t")

