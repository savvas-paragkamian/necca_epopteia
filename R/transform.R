## Script name: transform.R
##
## Authors: Savvas Paragkamian, Christina Kassara
##
## Date Created: 2025-07-01
##
## Purpose of script:
## Pure transform functions for the ETL pipeline: combines occurrence
## sources, assigns EEA grid cells, enriches records with elevation
## and border distance, applies species-specific distribution and
## population filters, and builds the final presence datasets.
## No file I/O; all paths are passed as arguments.
##
## Funded by: Natural Environment & Climate Change Agency (NECCA/ΟΦΥΠΕΚΑ)

library(sf)
library(terra)
library(units)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)
library(purrr)

build_national_report_distribution_grid_occurrences <- function(
  species_distribution,
  eea_grid_10km,
  species_taxonomy
) {
  species_dist_national_eea <- sf::st_join(
    eea_grid_10km,
    species_distribution,
    join = sf::st_within,
    left = FALSE,
    largest = FALSE
  )

  centroids_wgs <- sf::st_centroid(species_dist_national_eea) |>
    sf::st_transform(4326)
  coords <- sf::st_coordinates(centroids_wgs)

  species_dist_national_eea |>
    sf::st_drop_geometry() |>
    dplyr::mutate(
      decimalLongitude = coords[, 1],
      decimalLatitude  = coords[, 2]
    ) |>
    dplyr::distinct(submittedName, CELLCODE, decimalLongitude, decimalLatitude) |>
    dplyr::left_join(
      dplyr::select(species_taxonomy, verbatim_name, acceptedName),
      by = c("submittedName" = "verbatim_name")
    ) |>
    dplyr::rename(CELLCODE_eea_10km = CELLCODE, species = acceptedName) |>
    dplyr::filter(!is.na(species)) |>
    dplyr::mutate(
      collectionCode     = "DistrMap_2013_2018",
      basisOfRecord      = "ESTIMATED_CENTROID",
      recordNumber       = paste0("DistrMap_2013_2018_", sprintf("%04d", dplyr::row_number())),
      datasetName        = "GR_Art17_species_distribution.shp",
      DistrMap_2013_2018 = TRUE
    )
}

filter_e1x_ref_coords <- function(e1x_db_reference_occurrences) {
  dplyr::filter(e1x_db_reference_occurrences, !is.na(decimalLatitude))
}

build_e1x_ref_grid_occurrences <- function(e1x_db_reference_occurrences, eea_grid_10km) {
  cols <- c(occurrence_columns_to_keep(), "CELLCODE_eea_10km")

  records_no_coords <- e1x_db_reference_occurrences |>
    dplyr::filter(is.na(decimalLatitude)) |>
    dplyr::select(dplyr::all_of(cols))

  eea_10km_etrs89 <- sf::st_transform(eea_grid_10km, 3035)

  records_sf <- eea_10km_etrs89 |>
    dplyr::select(CELLCODE) |>
    dplyr::inner_join(records_no_coords, by = c("CELLCODE" = "CELLCODE_eea_10km"))

  centroids_wgs <- sf::st_centroid(records_sf) |>
    sf::st_transform(4326)
  coords <- sf::st_coordinates(centroids_wgs)

  records_sf |>
    sf::st_drop_geometry() |>
    dplyr::rename(CELLCODE_eea_10km = CELLCODE) |>
    dplyr::mutate(
      decimalLatitude  = coords[, 2],
      decimalLongitude = coords[, 1],
      basisOfRecord    = "MaterialCitation|ESTIMATED_CENTROID",
      collectionCode   = "E1X_DB_references_grid"
    )
}

build_p_apollo_grid_occurrences <- function(p_apollo_action_plan_occurrences) {
  centroids_wgs <- sf::st_centroid(p_apollo_action_plan_occurrences) |>
    sf::st_transform(4326)
  coords <- sf::st_coordinates(centroids_wgs)

  p_apollo_action_plan_occurrences |>
    sf::st_drop_geometry() |>
    dplyr::select(-dplyr::any_of(c("EOFORIGIN", "NOFORIGIN"))) |>
    dplyr::distinct() |>
    dplyr::mutate(
      decimalLongitude = coords[, 1],
      decimalLatitude  = coords[, 2],
      recordNumber     = paste0("action_plan_p.apollo_", sprintf("%02d", dplyr::row_number())),
      species          = "Parnassius apollo",
      individualCount  = NA_real_,
      submittedName    = "Parnassius apollo"
    )
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
  private_occurrences,
  p_apollo_action_plan_occurrences,
  e1x_ref_grid_occurrences
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
    private_occurrences,
    p_apollo_action_plan_occurrences,
    e1x_ref_grid_occurrences
  ) |>
    purrr::map(~ dplyr::select(.x, dplyr::all_of(cols))) |>
    dplyr::bind_rows()
}


normalize_taxonomy <- function(species_occurrences_invertebrates, species_taxonomy) {
  species_occurrences_invertebrates |>
    dplyr::left_join(
      dplyr::select(species_taxonomy, verbatim_name, acceptedName),
      by = c("submittedName" = "verbatim_name")
    ) |>
    dplyr::rename(species = acceptedName)
}

filter_art17_occurrences <- function(species_occurrences_normalized) {
  dplyr::filter(species_occurrences_normalized, !is.na(species))
}

# two levels of deduplication
# first is that the complete line is duplicated
# second is that the line is duplicated except from recordNumber. Meaning
# that the columns used here do not suffice for the uniquenness of the line.
deduplicate_art17_occurrences <- function(species_samples_art17) {
  species_samples_art17 |>
    dplyr::distinct() |>
    dplyr::group_by(dplyr::across(-recordNumber)) |> 
    dplyr::summarise(recordNumber = stringr::str_c(recordNumber, collapse = ", "),.groups = "keep") |>
    dplyr::ungroup()
}

filter_art17_coords <- function(species_samples_art17_dedup) {
  dplyr::filter(species_samples_art17_dedup, !is.na(decimalLatitude))
}

assign_eea_grid_10km <- function(species_samples_art17_coords, eea_grid_10km) {
  eea_10km_etrs89 <- sf::st_transform(eea_grid_10km, 3035)

  species_samples_art17_coords |>
    sf::st_as_sf(
      coords = c("decimalLongitude", "decimalLatitude"),
      remove = FALSE,
      crs    = "WGS84"
    ) |>
    sf::st_transform(3035) |>
    sf::st_join(eea_10km_etrs89) |>
    sf::st_drop_geometry() |>
    dplyr::rename(CELLCODE_eea_10km = CELLCODE) |>
    dplyr::distinct(submittedName, decimalLatitude, decimalLongitude, collectionCode,
                    datasetName, recordNumber, basisOfRecord, individualCount,
                    CELLCODE_eea_10km, species)
}

find_distrmap_cells_without_e1x <- function(
  species_samples_eea,
  national_report_distribution_grid
) {
  e1x_db_cells <- species_samples_eea |>
    dplyr::filter(stringr::str_detect(collectionCode, "^E1X_DB")) |>
    dplyr::distinct(species, CELLCODE_eea_10km)

  national_report_distribution_grid |>
    dplyr::distinct(species, CELLCODE_eea_10km) |>
    dplyr::anti_join(e1x_db_cells, by = c("species", "CELLCODE_eea_10km"))
}

find_national_report_orphan_cells <- function(
  species_samples_eea,
  national_report_distribution_grid
) {
  eea_cells_no_gbif <- species_samples_eea |>
    dplyr::filter(collectionCode != "GBIF") |>
    dplyr::distinct(species, CELLCODE_eea_10km)

  national_report_distribution_grid |>
    dplyr::anti_join(eea_cells_no_gbif, by = c("species", "CELLCODE_eea_10km")) |>
    dplyr::mutate(individualCount = NA_real_)
}

filter_e1x_ref_grid <- function(species_samples_eea) {
  dplyr::filter(species_samples_eea, collectionCode != "E1X_DB_references_grid")
}

build_presence_minimum <- function(
  species_samples_eea,
  national_report_distribution_grid,
  national_report_orphan_cells
) {
  distrmap_flag <- national_report_distribution_grid |>
    dplyr::distinct(species, CELLCODE_eea_10km, DistrMap_2013_2018)

  species_samples_flagged <- species_samples_eea |>
    dplyr::left_join(distrmap_flag, by = c("species", "CELLCODE_eea_10km")) |>
    dplyr::mutate(
      DistrMap_2013_2018 = dplyr::if_else(is.na(DistrMap_2013_2018), FALSE, DistrMap_2013_2018)
    )

  dplyr::bind_rows(species_samples_flagged, national_report_orphan_cells)
}

enrich_with_elevation <- function(species_samples_presence_minimum, eu_dem_path) {
  eu_dem <- terra::rast(eu_dem_path)

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
      collectionCode == "GBIF",
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      collectionCode == "E1X_DB_references" & !DistrMap_2013_2018,
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
      species == "Parnassius apollo" & eudem_dem_3035_europe <= 450 & collectionCode != "Action Plan 2019",
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      species == "Parnassius apollo" & collectionCode %in% c("E1X_DB_references", "DistrMap_2013_2018"),
      FALSE, includeDistribution)) |>
    dplyr::mutate(includeDistribution = dplyr::if_else(
      species == "Rosalia alpina" & minimumDistanceFromBorders > 500 & collectionCode == "NECCA_redlist",
      TRUE, includeDistribution))
}

apply_population_filters <- function(species_samples_presence_dist, eea_grid_1km) {
  eea_1km_etrs89 <- sf::st_transform(eea_grid_1km, 3035) |>
    dplyr::select(CELLCODE)

  species_samples_presence_dist |>
    sf::st_join(eea_1km_etrs89) |>
    dplyr::rename(CELLCODE_eea_1km = CELLCODE) |>
    dplyr::mutate(includePopulation = dplyr::if_else(
      collectionCode %in% c("E1X_DB_references", "DistrMap_2013_2018"),
      FALSE, includeDistribution)) |>
    dplyr::mutate(includePopulation = dplyr::if_else(
      species == "Vertigo angustior", FALSE, includePopulation))
}

build_presence_final <- function(species_samples_presence_pop) {
  species_samples_presence_pop |>
    sf::st_drop_geometry() |>
    dplyr::filter(collectionCode != "Invertebrates_records_private")
}

build_presence_final_private <- function(species_samples_presence_pop) {
  sf::st_drop_geometry(species_samples_presence_pop)
}

build_presence_final_sf <- function(species_samples_presence_pop) {
  dplyr::filter(species_samples_presence_pop,
                collectionCode != "Invertebrates_records_private")
}

compute_species_range <- function(species_samples_presence_final, eea_grid_10km) {
  eea_10km_etrs89 <- sf::st_transform(eea_grid_10km, 3035)

  presence <- dplyr::filter(species_samples_presence_final,
                             includeDistribution == TRUE)

  species_list <- sort(unique(presence$species))

  range_list <- lapply(species_list, function(sp) {
    sp_cells <- presence |>
      dplyr::filter(species == sp) |>
      dplyr::distinct(CELLCODE_eea_10km) |>
      dplyr::rename(CELLCODE = CELLCODE_eea_10km)

    grids    <- dplyr::filter(eea_10km_etrs89, CELLCODE %in% sp_cells$CELLCODE)
    expanded <- expand_range_with_gap_distance(
      distribution   = grids,
      full_grid      = eea_10km_etrs89,
      gap_distance_m = 40000
    )
    dplyr::mutate(expanded, species = sp)
  })

  dplyr::bind_rows(range_list)
}

enrich_with_natura2000 <- function(species_samples_presence_final_sf, natura2000) {
  natura2000_etrs89 <- sf::st_transform(natura2000, 3035)
  sf::st_join(
    species_samples_presence_final_sf,
    dplyr::select(natura2000_etrs89, SITECODE, SITETYPE)
  )
}

build_distributions_summary <- function(
  species_samples_presence_final,
  eea_grid_10km,
  greece_regions
) {
  eea_10km_etrs89 <- sf::st_transform(eea_grid_10km, 3035)
  greece_etrs89   <- sf::st_transform(greece_regions, 3035)

  # Land area per cell: intersection clips coastal/border cells to actual territory
  cell_area_land <- sf::st_intersection(eea_10km_etrs89, greece_etrs89) |>
    dplyr::mutate(AreaKm2 = as.numeric(units::set_units(sf::st_area(geometry), "km^2"))) |>
    sf::st_drop_geometry() |>
    dplyr::group_by(CELLCODE) |>
    dplyr::summarise(AreaKm2 = sum(AreaKm2), .groups = "drop")

  species_samples_presence_final |>
    dplyr::filter(includeDistribution) |>
    dplyr::distinct(species, CELLCODE_eea_10km) |>
    dplyr::rename(CELLCODE = CELLCODE_eea_10km) |>
    dplyr::left_join(cell_area_land, by = "CELLCODE") |>
    dplyr::group_by(species) |>
    dplyr::summarise(
      CellCount = dplyr::n(),
      AreaKm2   = round(sum(AreaKm2, na.rm = TRUE), 2),
      .groups   = "drop"
    )
}

build_populations_summary <- function(species_samples_presence_final) {
  populations_with_refs <- species_samples_presence_final |>
    dplyr::filter(includeDistribution) |>
    dplyr::distinct(species, CELLCODE_eea_1km) |>
    dplyr::group_by(species) |>
    dplyr::summarise(n_1km = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(method = "with_refs")

  populations_no_refs <- species_samples_presence_final |>
    dplyr::filter(includePopulation, collectionCode != "E1X_DB_references") |>
    dplyr::distinct(species, CELLCODE_eea_1km) |>
    dplyr::group_by(species) |>
    dplyr::summarise(n_1km = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(method = "no_refs")

  dplyr::bind_rows(populations_with_refs, populations_no_refs) |>
    tidyr::pivot_wider(names_from = method, values_from = n_1km)
}

build_natura_populations_summary <- function(species_samples_presence_natura) {
  species_samples_presence_natura |>
    sf::st_drop_geometry() |>
    dplyr::filter(includePopulation, !is.na(SITECODE)) |>
    dplyr::distinct(species, SITECODE, SITETYPE, CELLCODE_eea_1km) |>
    dplyr::group_by(species, SITECODE, SITETYPE) |>
    dplyr::summarise(n_1km = dplyr::n(), .groups = "drop")
}
