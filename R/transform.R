library(sf)
library(terra)
library(units)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(purrr)

build_national_report_distribution_grid_occurrences <- function(
  species_distribution,
  eea_grid_10km
) {
  species_dist_national_eea <- sf::st_join(
    eea_grid_10km,
    species_distribution,
    join = sf::st_within,
    left = FALSE,
    largest = FALSE
  )

  centroids <- sf::st_centroid(species_dist_national_eea)
  coords    <- sf::st_coordinates(centroids)

  species_dist_national_eea <- cbind(species_dist_national_eea,
                                     decimalLongitude = coords[, 1],
                                     decimalLatitude  = coords[, 2]) |>
    dplyr::distinct(submittedName, species, CELLCODE, decimalLongitude, decimalLatitude) |>
    dplyr::rename(CELLCODE_eea_10km = CELLCODE) |>
    dplyr::mutate(
      datasetName    = "DistrMap_2013_2018",
      basisOfRecord  = "ESTIMATED_CENTROID",
      recordNumber   = paste0("DistrMap_2013_2018_", sprintf("%04d", dplyr::row_number())),
      collectionCode = "GR_Art17_species_distribution.shp",
      DistrMap_2013_2018 = TRUE
    )

  grid_laea <- species_dist_national_eea |>
    sf::st_as_sf(
      coords = c("decimalLongitude", "decimalLatitude"),
      remove = TRUE,
      crs    = sf::st_crs(species_distribution)
    )

  grid_wgs <- sf::st_transform(grid_laea, 4326)
  coords_wgs <- sf::st_coordinates(grid_wgs)

  cbind(grid_wgs,
        decimalLongitude = coords_wgs[, 1],
        decimalLatitude  = coords_wgs[, 2])
}

combine_all_occurrences <- function(
  gbif_occurrences,
  e1x_mdpp_occurrences,
  e1x_db_sampling_occurrences,
  e1x_db_reference_occurrences,
  e2x_occurrences,
  necca_redlist_points_occurrences,
  e2x_ref_unio_crassus,
  e2x_ref_stenobothrus_eurasius,
  private_occurrences
) {
  cols <- occurrence_columns_to_keep()
  list(
    gbif_occurrences,
    e1x_mdpp_occurrences,
    e1x_db_sampling_occurrences,
    e1x_db_reference_occurrences,
    e2x_occurrences,
    necca_redlist_points_occurrences,
    e2x_ref_unio_crassus,
    e2x_ref_stenobothrus_eurasius,
    private_occurrences
  ) |>
    purrr::map(~ dplyr::select(.x, dplyr::all_of(cols))) |>
    dplyr::bind_rows()
}

filter_art17_occurrences <- function(
  species_occurrences_invertebrates,
  species_taxonomy
) {
  species_names <- unique(species_taxonomy$verbatim_name)
  species_occurrences_invertebrates |>
    dplyr::filter(submittedName %in% species_names) |>
    dplyr::left_join(species_taxonomy, by = c("submittedName" = "verbatim_name"))
}

assign_eea_grid_10km <- function(species_samples_art17, eea_grid_10km) {
  eea_10km_etrs89 <- sf::st_transform(eea_grid_10km, 3035)

  species_samples_art17 |>
    dplyr::filter(!is.na(decimalLatitude)) |>
    dplyr::distinct(submittedName, decimalLatitude, decimalLongitude, datasetName,
                    collectionCode, recordNumber, basisOfRecord, individualCount, species) |>
    dplyr::group_by(dplyr::across(-recordNumber)) |>
    dplyr::summarise(recordNumber = stringr::str_c(recordNumber, collapse = ", "),
                     .groups = "keep") |>
    dplyr::mutate(species = dplyr::case_when(
      species == "Hirudo medicinalis"    ~ "Hirudo verbana",
      species == "Unio vicarius"         ~ "Unio crassus",
      species == "Unio ionicus"          ~ "Unio crassus",
      species == "Unio bruguierianus"    ~ "Unio crassus",
      species == "Unio desectus"         ~ "Unio crassus",
      species == "Polyommatus eros"      ~ "Polyommatus eroides",
      species == "Paracossulus thrips"   ~ "Catopta thrips",
      species == "Osmoderma eremita"     ~ "Osmoderma lassallei",
      .default = species
    )) |>
    sf::st_as_sf(
      coords = c("decimalLongitude", "decimalLatitude"),
      remove = FALSE,
      crs    = "WGS84"
    ) |>
    sf::st_transform(3035) |>
    sf::st_join(eea_10km_etrs89) |>
    sf::st_drop_geometry() |>
    dplyr::rename(CELLCODE_eea_10km = CELLCODE) |>
    dplyr::distinct(submittedName, decimalLatitude, decimalLongitude, datasetName,
                    collectionCode, recordNumber, basisOfRecord, individualCount,
                    CELLCODE_eea_10km, species)
}

build_presence_minimum <- function(
  species_samples_eea,
  p_apollo_action_plan_occurrences,
  national_report_distribution_grid
) {
  # Combine samples with Apollo action plan (both already have CELLCODE_eea_10km)
  species_samples_combined <- dplyr::bind_rows(
    species_samples_eea,
    p_apollo_action_plan_occurrences
  )

  # Which cells appear in the 2013-2018 national report
  species_dist_national_minimum <- national_report_distribution_grid |>
    sf::st_drop_geometry() |>
    dplyr::distinct(species, CELLCODE_eea_10km, DistrMap_2013_2018)

  # Attach DistrMap flag and composite key for orphan detection
  species_samples_combined_dist <- species_samples_combined |>
    dplyr::left_join(species_dist_national_minimum,
                     by = c("species", "CELLCODE_eea_10km")) |>
    dplyr::mutate(
      DistrMap_2013_2018 = dplyr::if_else(is.na(DistrMap_2013_2018), FALSE, DistrMap_2013_2018),
      composite_key      = paste0(species, "_", CELLCODE_eea_10km)
    )

  # National report as flat dataset (geometry and DistrMap flag dropped)
  national_report_dataset <- national_report_distribution_grid |>
    sf::st_drop_geometry() |>
    dplyr::select(-DistrMap_2013_2018) |>
    dplyr::mutate(composite_key = paste0(species, "_", CELLCODE_eea_10km))

  # Orphan cells: in national report but absent from non-GBIF point data
  samples_no_gbif_keys <- species_samples_combined_dist |>
    dplyr::filter(datasetName != "GBIF") |>
    dplyr::pull(composite_key) |>
    unique()

  national_orphans <- national_report_dataset |>
    dplyr::filter(!composite_key %in% samples_no_gbif_keys) |>
    dplyr::mutate(DistrMap_2013_2018 = TRUE, individualCount = NA)

  dplyr::bind_rows(species_samples_combined_dist, national_orphans) |>
    dplyr::select(-composite_key)
}

enrich_with_elevation <- function(species_samples_presence_minimum, eu_dem) {
  presence_sf <- species_samples_presence_minimum |>
    sf::st_as_sf(
      coords = c("decimalLongitude", "decimalLatitude"),
      remove = FALSE,
      crs    = "WGS84"
    )

  vect_points <- terra::vect(presence_sf)
  vals        <- terra::extract(eu_dem, vect_points)[, -1, drop = FALSE]
  names(vals) <- "eudem_dem_3035_europe"

  cbind(presence_sf, vals)
}

enrich_with_border_distance <- function(species_samples_presence_sf, greece_regions) {
  presence_etrs89 <- sf::st_transform(species_samples_presence_sf, 3035)
  greece_etrs89   <- sf::st_transform(greece_regions, 3035)
  dist_matrix     <- sf::st_distance(presence_etrs89, greece_etrs89)
  presence_etrs89$minimumDistanceFromBorders <- apply(dist_matrix, 1, min)
  presence_etrs89
}

apply_distribution_filters <- function(species_samples_presence_etrs89) {
  species_samples_presence_etrs89 |>
    dplyr::mutate(
      species = dplyr::if_else(
        species == "Zerynthia polyxena" & decimalLatitude < 36,
        "Zerynthia cretica", species
      )
    ) |>
    dplyr::mutate(includeDistribution = TRUE) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      datasetName == "GBIF",
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      datasetName == "E1X_DB_references" & !DistrMap_2013_2018,
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      species == "Zerynthia cretica",
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      species == "Vertigo moulinsiana",
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      minimumDistanceFromBorders > 500 & basisOfRecord != "ESTIMATED_CENTROID",
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      species == "Stenobothrus eurasius" & CELLCODE_eea_10km == "10kmE535N184",
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      species == "Probaticus subrugosus" & CELLCODE_eea_10km == "10kmE524N175",
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      species == "Pseudophilotes bavius" & CELLCODE_eea_10km == "10kmE535N176",
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      species == "Parnassius apollo" & eudem_dem_3035_europe <= 450 & datasetName != "Action Plan 2019",
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      species == "Parnassius apollo" & datasetName %in% c("E1X_DB_references", "DistrMap_2013_2018"),
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      species == "Rosalia alpina" & minimumDistanceFromBorders > 500 & datasetName == "NECCA_redlist",
      TRUE, includeDistribution))
}

apply_population_filters <- function(species_samples_presence_dist, eea_grid_1km) {
  eea_1km_etrs89 <- sf::st_transform(eea_grid_1km, 3035) |>
    dplyr::select(CELLCODE)

  species_samples_presence_dist |>
    sf::st_join(eea_1km_etrs89) |>
    dplyr::rename(CELLCODE_eea_1km = CELLCODE) |>
    dplyr::mutate(includePopulation = dplyr::if_else(
      datasetName %in% c("E1X_DB_references", "DistrMap_2013_2018"),
      FALSE, includeDistribution)) |>
    dplyr::mutate(includePopulation = dplyr::if_else(
      species == "Vertigo angustior", FALSE, includePopulation))
}

build_presence_final <- function(species_samples_presence_pop) {
  species_samples_presence_pop |>
    sf::st_drop_geometry() |>
    dplyr::filter(datasetName != "Invertebrates_records_private")
}

build_presence_final_private <- function(species_samples_presence_pop) {
  sf::st_drop_geometry(species_samples_presence_pop)
}
