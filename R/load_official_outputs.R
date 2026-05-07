## Script name: load_official_outputs.R
##
## Authors: Savvas Paragkamian, Christina Kassara
##
## Date Created: 2025-07-01
##
## Purpose of script:
## TSV writers for the official Article 17 reporting files:
## combined occurrences, Art. 17 filtered records, final presence
## dataset, per-species distributions (10 km grid, observed +
## range), and per-species population counts (1 km grid).
## Each function returns the written file path for targets
## format = "file" tracking.
##
## Funded by: Natural Environment & Climate Change Agency (NECCA/ΟΦΥΠΕΚΑ)

library(dplyr)
library(readr)
library(sf)

save_occurrences_tsv <- function(species_occurrences_invertebrates, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_delim(species_occurrences_invertebrates, path, delim = "\t")
  path
}

save_art17_tsv <- function(species_samples_art17, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_delim(species_samples_art17, path, delim = "\t")
  path
}

save_presence_final_tsv <- function(species_samples_presence_final, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_delim(species_samples_presence_final, path, delim = "\t")
  path
}

save_distributions_tsv <- function(
  species_samples_presence_final,
  species_range,
  path
) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  # Observation-based: unique (species, 10km cell) where includeDistribution is TRUE
  obs_distribution <- species_samples_presence_final |>
    dplyr::filter(includeDistribution) |>
    dplyr::group_by(species, CELLCODE_eea_10km) |>
    dplyr::summarise(
      in_observation = TRUE,
      n_records      = dplyr::n(),
      n_datasets     = dplyr::n_distinct(datasetName),
      datasets       = paste(sort(unique(datasetName)), collapse = ";"),
      .groups        = "drop"
    )

  # Range cells from expand_range_with_gap_distance (includes inferred gap-fill cells)
  range_cells <- species_range |>
    sf::st_drop_geometry() |>
    dplyr::distinct(species, CELLCODE, cell_origin) |>
    dplyr::rename(CELLCODE_eea_10km = CELLCODE) |>
    dplyr::mutate(in_range = TRUE)

  # Full distribution: observed cells union inferred range cells
  distributions <- dplyr::full_join(
    obs_distribution,
    range_cells,
    by = c("species", "CELLCODE_eea_10km")
  ) |>
    dplyr::mutate(
      in_observation = dplyr::if_else(is.na(in_observation), FALSE, in_observation),
      in_range       = dplyr::if_else(is.na(in_range),       FALSE, in_range)
    )

  readr::write_delim(distributions, path, delim = "\t")
  path
}

save_populations_tsv <- function(species_samples_presence_final, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  populations <- species_samples_presence_final |>
    dplyr::filter(includePopulation) |>
    dplyr::group_by(species, CELLCODE_eea_10km, CELLCODE_eea_1km) |>
    dplyr::summarise(
      n_records         = dplyr::n(),
      n_datasets        = dplyr::n_distinct(datasetName),
      total_individuals = sum(as.numeric(individualCount), na.rm = TRUE),
      datasets          = paste(sort(unique(datasetName)), collapse = ";"),
      .groups           = "drop"
    )

  readr::write_delim(populations, path, delim = "\t")
  path
}
