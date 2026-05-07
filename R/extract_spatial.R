## Script name: extract_spatial.R
##
## Authors: Savvas Paragkamian
##
## Date Created: 2025-02-01
##
## Purpose of script:
## Loads all spatial reference layers (Greece boundaries, EEA 1 km
## and 10 km grids, Natura 2000 sites) into a named list for use
## throughout the pipeline. The EU DEM raster is excluded here and
## read on demand in transform.R (terra external pointer constraint).
##
## Funded by: Natural Environment & Climate Change Agency (NECCA/ΟΦΥΠΕΚΑ)

library(sf)
library(terra)

read_spatial_layers <- function(paths) {
  list(
    greece_regions = sf::st_read(paths$greece_regions, quiet = TRUE),
    eea_grid_1km   = sf::st_read(paths$eea_grid_1km, quiet = TRUE),
    eea_grid_10km  = sf::st_read(paths$eea_grid_10km, quiet = TRUE),
    natura2000     = sf::st_read(paths$natura2000, quiet = TRUE)
  )
}
