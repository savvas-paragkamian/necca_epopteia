## Script name: extract_occurrences.R
##
## Authors: Savvas Paragkamian, Christina Kassara
##
## Date Created: 2024-11-01
##
## Purpose of script:
## One reader function per occurrence data source (GBIF, E1X MDPP,
## E1X DB sampling, E1X DB references, E2X DB, private records,
## Unio crassus supplement, Stenobothrus eurasius supplement,
## NECCA Redlist, national report distribution, Parnassius apollo
## action plan). All outputs use Darwin Core column names.
##
## Funded by: Natural Environment & Climate Change Agency (NECCA/ΟΦΥΠΕΚΑ)

library(readxl)
library(taxize)
library(rgbif)
library(dplyr)
library(readr)


read_species_taxonomy <- function(path) {
  readr::read_delim(path, delim = "\t", show_col_types = FALSE)
}

get_species_names_combined <- function(species_taxonomy) {
  as.character(species_taxonomy$verbatim_name)
}

occurrence_columns_to_keep <- function() {
  c(
    "submittedName",
    "decimalLatitude",
    "decimalLongitude",
    "datasetName",
    "recordNumber",
    "collectionCode",
    "basisOfRecord",
    "individualCount"
  )
}

read_gbif_occurrences <- function(path, greece_regions) {
  gbif_species_occ <- readr::read_delim(path, delim = "\t", show_col_types = FALSE) |>
    dplyr::mutate(datasetName = "GBIF") |>
    dplyr::mutate(
      species = ifelse(
        !is.na(verbatimScientificName) &
          verbatimScientificName == "Panaxia quadripunctaria",
        "Euplagia quadripunctaria",
        as.character(species)
      )
    ) |>
    dplyr::rename(submittedName = species) |>
    dplyr::filter(coordinateUncertaintyInMeters < 1000)

  gbif_species_occ_sf <- gbif_species_occ |>
    sf::st_as_sf(
      coords = c("decimalLongitude", "decimalLatitude"),
      remove = FALSE,
      crs = "WGS84"
    )

  within_mat <- sf::st_intersects(gbif_species_occ_sf, greece_regions, sparse = FALSE)

  gbif_species_occ_sf <- gbif_species_occ_sf[rowSums(within_mat) > 0, ]

  gbif_species_occ_gr <- gbif_species_occ_sf |>
    sf::st_drop_geometry() |>
    dplyr::mutate(recordNumber = occurrenceID) |>
    dplyr::mutate(collectionCode = basename(path))

  gbif_species_occ_gr
}

read_e1x_mdpp_occurrences <- function(path) {
  samples_data <- readxl::read_xlsx(
    path,
    sheet = "Δείγματα Ασπόνδυλων",
    col_names = TRUE
  ) |>
    dplyr::slice(-1) |>
    dplyr::mutate(
      decimalLatitude = as.numeric(`Γεωγραφικό Πλάτος (WGS84) Αρχη`),
      decimalLongitude = as.numeric(`Γεωγραφικό Μήκος (WGS84) Αρχή`)
    ) |>
    dplyr::filter(!is.na(decimalLongitude)) |>
    dplyr::mutate(datasetName = "E1X_MDPP_2014_2024") |>
    dplyr::mutate(basisOfRecord = "MATERIAL_SAMPLE")

  species_data <- readxl::read_xlsx(
    path,
    sheet = "Είδη",
    col_names = TRUE
  ) |>
    dplyr::slice(-1)

  species_data |>
    dplyr::mutate(
      submittedName = dplyr::if_else(`Όνομα είδους` == "Άλλο", `Άλλο είδος`, `Όνομα είδους`)
    ) |>
    dplyr::mutate(
      art17_92_43_EEC = dplyr::if_else(`Όνομα είδους` != "Άλλο", TRUE, FALSE)
    ) |>
    dplyr::mutate(individualCount = as.numeric(`Αριθμός ατόμων είδους`)) |>
    dplyr::mutate(organismQuantity = as.numeric(`Κατηγορία Σχετικής αφθονίας είδους`)) |>
    dplyr::filter(organismQuantity != 0 | is.na(organismQuantity)) |>
    dplyr::left_join(samples_data, by = c("Sam_ID" = "Sam_ID")) |>
    dplyr::mutate(recordNumber = Obs_ID) |>
    dplyr::mutate(collectionCode = basename(path))
}

read_e1x_db_reference_occurrences <- function(path) {
  ref_samples_data <- readxl::read_xlsx(
    path,
    sheet = "Εξάπλωση ειδών και τ.ο.",
    col_names = TRUE
  ) |>
    dplyr::slice(-1) |>
    dplyr::filter(!is.na(`Κωδικός Αναφοράς`))

  refs_data <- readxl::read_xlsx(
    path,
    sheet = "Βιβλιογραφία",
    col_names = TRUE
  ) |>
    dplyr::slice(-1) |>
    dplyr::filter(!is.na(`Κωδικός Αναφοράς`))

  ref_samples_data$decimalLatitude <- as.numeric(ref_samples_data$`Γεωγραφικό Πλάτος (WGS84)`)
  ref_samples_data$decimalLongitude <- as.numeric(ref_samples_data$`Γεωγραφικό Μήκος (WGS84)`)
  ref_samples_data$species <- ref_samples_data$`Ονομασία είδους`

  ref_samples_data |>
    dplyr::left_join(refs_data, by = c("Κωδικός Αναφοράς" = "Κωδικός Αναφοράς")) |>
    dplyr::mutate(datasetName = "E1X_DB_references") |>
    dplyr::mutate(basisOfRecord = "MaterialCitation") |>
    dplyr::mutate(submittedName = `Ονομασία είδους`) |>
    dplyr::mutate(individualCount = as.numeric(`Πλήθος ατόμων`)) |>
    dplyr::mutate(recordNumber = SpRef_ID) |>
    dplyr::mutate(collectionCode = basename(path))
}

read_e1x_db_sampling_occurrences <- function(path) {
  samples_data <- readxl::read_xlsx(
    path,
    sheet = "Δείγματα Ασπόνδυλων",
    col_names = TRUE
  ) |>
    dplyr::slice(-1)

  samples_data$decimalLatitude <- as.numeric(samples_data$`Γεωγραφικό Πλάτος (WGS84) Αρχη`)
  samples_data$decimalLongitude <- as.numeric(samples_data$`Γεωγραφικό Μήκος (WGS84) Αρχή`)

  species_data <- readxl::read_xlsx(
    path,
    sheet = "Είδη",
    col_names = TRUE
  ) |>
    dplyr::slice(-1)

  all_data <- species_data |>
    dplyr::filter(!is.na(Sam_ID)) |>
    dplyr::left_join(samples_data, by = c("Sam_ID" = "Sam_ID")) |>
    dplyr::mutate(
      submittedName = dplyr::if_else(`Όνομα είδους` == "Άλλο", `Άλλο είδος`, `Όνομα είδους`)
    ) |>
    dplyr::mutate(
      art17_92_43_EEC = dplyr::if_else(`Όνομα είδους` != "Άλλο", TRUE, FALSE)
    ) |>
    dplyr::mutate(individualCount = as.numeric(`Αριθμός ατόμων είδους`)) |>
    dplyr::mutate(datasetName = "E1X_DB") |>
    dplyr::mutate(basisOfRecord = "MATERIAL_SAMPLE") |>
    dplyr::mutate(recordNumber = Obs_ID) |>
    dplyr::mutate(collectionCode = basename(path))

  all_data |>
    dplyr::filter(individualCount > 0)
}

read_e2x_occurrences <- function(path) {
  readr::read_delim(path, delim = "\t", show_col_types = FALSE) |>
    dplyr::mutate(individualCount = NA) |>
    dplyr::mutate(recordNumber = paste0("E2X_DB_", sprintf("%02d", dplyr::row_number()))) |>
    dplyr::mutate(collectionCode = basename(path))
}

read_private_occurrences <- function(path) {
  readr::read_delim(path, delim = ",", show_col_types = FALSE) |>
    dplyr::mutate(
      decimalLongitude = Longitude,
      decimalLatitude = Latitude,
      submittedName = Species,
      individualCount = as.numeric(Individuals)
    ) |>
    dplyr::mutate(datasetName = "Invertebrates_records_private") |>
    dplyr::mutate(basisOfRecord = "MATERIAL_SAMPLE") |>
    dplyr::mutate(recordNumber = as.character(ID)) |>
    dplyr::bind_rows(
      tibble::tibble(
        submittedName = "Cerambyx cerdo",
        decimalLatitude = 38.077228,
        decimalLongitude = 24.380490,
        datasetName = "Invertebrates_records_private",
        recordNumber = "157",
        basisOfRecord = "MATERIAL_SAMPLE"
      )
    ) |>
    dplyr::mutate(collectionCode = basename(path))
}

read_unio_crassus_occurrences <- function(path) {
  readxl::read_xlsx(path, col_names = TRUE) |>
    dplyr::rename(
      decimalLongitude = LONGITUDE,
      decimalLatitude = LATITUDE,
      submittedName = SPECIES
    ) |>
    dplyr::mutate(
      submittedName = gsub("U.", "Unio", submittedName),
      datasetName = "E2X_ref",
      collectionCode = basename(path),
      recordNumber = paste0("Lopes-Lima_2024_", sprintf("%02d", dplyr::row_number())),
      basisOfRecord = "MaterialCitation"
    ) |>
    dplyr::mutate(individualCount = NA) |>
    dplyr::filter(COUNTRY == "Greece") |>
    dplyr::distinct(
      datasetName,
      basisOfRecord,
      submittedName,
      collectionCode,
      recordNumber,
      individualCount,
      decimalLatitude,
      decimalLongitude
    )
}

read_stenobothrus_occurrences <- function(path) {
  readr::read_delim(path, delim = "\t", show_col_types = FALSE)
}

read_necca_redlist_points_occurrences <- function(path) {
  x <- sf::st_read(path, quiet = TRUE) |>
    sf::st_cast("POINT") 
  coords <- sf::st_coordinates(x)
  
  x <- x |>
    dplyr::mutate(
      decimalLongitude = coords[, 1],
      decimalLatitude = coords[, 2]
    ) |>
    sf::st_drop_geometry() |>
    dplyr::rename(submittedName = sci_name) |>
    dplyr::mutate(datasetName = "NECCA_redlist") |>
    dplyr::mutate(individualCount = NA) |>
    dplyr::mutate(basisOfRecord = "MATERIAL_SAMPLE") |>
    dplyr::mutate(recordNumber = paste0("NECCA_redlist_", sprintf("%03d", dplyr::row_number()))) |>
    dplyr::mutate(collectionCode = basename(path))
}

read_national_report_distribution_occurrences <- function(
  distribution_path,
  sensitive_distribution_path,
  species_taxonomy
) {
  species_names_combined <- as.character(species_taxonomy$verbatim_name)

  species_dist_national_rep <- sf::st_read(
    distribution_path,
    quiet = TRUE
  )

  species_dist_national_rep_sens <- sf::st_read(
    sensitive_distribution_path,
    quiet = TRUE
  )

  dplyr::bind_rows(
    species_dist_national_rep,
    species_dist_national_rep_sens
  ) |>
    sf::st_make_valid() |>
    dplyr::mutate(
      submittedName = gsub("\\* ?", "", iname),
      submittedName = gsub(" Complex", "", submittedName),
      submittedName = gsub(
        "Osmoderma eremita",
        "Osmoderma lassallei",
        submittedName
      )
    ) |>
    dplyr::filter(submittedName %in% species_names_combined) |>
    dplyr::left_join(
      species_taxonomy,
      by = c("submittedName" = "verbatim_name")
    ) |>
    dplyr::mutate(
      species = gsub("Osmoderma eremita", "Osmoderma lassallei", species),
      species = gsub("Hirudo medicinalis", "Hirudo verbana", species),
      species = gsub("Unio vicarius", "Unio crassus", species),
      species = gsub("Unio ionicus", "Unio crassus", species),
      species = gsub("Unio bruguierianus", "Unio crassus", species),
      species = gsub("Unio desectus", "Unio crassus", species),
      species = gsub("Polyommatus eros", "Polyommatus eroides", species),
      species = gsub("Paracossulus thrips", "Catopta thrips", species),
      datasetName = "National_report_2013_2018",
      basisOfRecord = "ESTIMATED_DISTRIBUTION",
      collectionCode = paste(
        basename(distribution_path),
        basename(sensitive_distribution_path),
        sep = ";"
      )
    )
}

read_p_apollo_action_plan_occurrences <- function(path, eea_grid_10km) {

  p_apollo_dist <- sf::st_read(path, quiet = TRUE) |>
    dplyr::mutate(
      datasetName = "Action Plan 2019",
      basisOfRecord = "ESTIMATED_CENTROID"
    )

  p_apollo_centroids <- sf::st_centroid(p_apollo_dist)
  p_apollo_coords <- sf::st_coordinates(p_apollo_centroids)

  p_apollo_eea_coords <- cbind(
    p_apollo_dist,
    decimalLongitude = p_apollo_coords[, 1],
    decimalLatitude = p_apollo_coords[, 2]
  )

  p_apollo_coords_laea <- p_apollo_eea_coords |>
    sf::st_drop_geometry() |>
    sf::st_as_sf(
      coords = c("decimalLongitude", "decimalLatitude"),
      remove = TRUE,
      crs = sf::st_crs(eea_grid_10km)
    )

  p_apollo_coords_wgs <- sf::st_transform(
    p_apollo_coords_laea,
    4326
  )

  p_apollo_coords_d <- sf::st_coordinates(p_apollo_coords_wgs)

  p_apollo_eea_coords_wgs <- cbind(
    p_apollo_coords_wgs,
    decimalLongitude = p_apollo_coords_d[, 1],
    decimalLatitude = p_apollo_coords_d[, 2]
  )

  p_apollo_eea_coords_wgs |>
    sf::st_drop_geometry() |>
    dplyr::select(
      -dplyr::any_of(c("EOFORIGIN", "NOFORIGIN"))
    ) |>
    dplyr::distinct() |>
    dplyr::rename(
      CELLCODE_eea_10km = CELLCODE
    ) |>
    dplyr::mutate(
      recordNumber = paste0(
        "action_plan_p.apollo_",
        sprintf("%02d", dplyr::row_number())
      ),
      collectionCode = basename(path),
      species = "Parnassius apollo",
      individualCount = NA,
      submittedName = "Parnassius apollo"
    )
}
