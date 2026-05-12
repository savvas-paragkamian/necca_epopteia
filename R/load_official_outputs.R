## Script name: load_official_outputs.R
##
## Authors: Savvas Paragkamian, Christina Kassara
##
## Date Created: 2025-07-01
##
## Purpose of script:
## Thin TSV writers for the official Article 17 reporting files.
## Each function receives a ready-to-write data frame, writes it,
## and returns the file path for targets format = "file" tracking.
## All aggregation logic lives in transform.R (build_*_summary).
##
## Funded by: Natural Environment & Climate Change Agency (NECCA/ΟΦΥΠΕΚΑ)

library(readr)

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

save_distributions_tsv <- function(distributions_summary, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_delim(distributions_summary, path, delim = "\t")
  path
}

save_populations_tsv <- function(populations_summary, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_delim(populations_summary, path, delim = "\t")
  path
}

save_natura_populations_tsv <- function(natura_populations_summary, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_delim(natura_populations_summary, path, delim = "\t")
  path
}
