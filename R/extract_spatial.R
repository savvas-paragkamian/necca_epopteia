library(sf)
library(terra)

read_spatial_layers <- function(paths) {
  list(
    greece_regions = sf::st_read(paths$greece_regions, quiet = TRUE),
    eea_grid_1km   = sf::st_read(paths$eea_grid_1km, quiet = TRUE),
    eea_grid_10km  = sf::st_read(paths$eea_grid_10km, quiet = TRUE),
    natura2000     = sf::st_read(paths$natura2000, quiet = TRUE),
    eu_dem         = terra::rast(paths$eu_dem)
  )
}
