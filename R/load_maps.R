## Script name: load_maps.R
##
## Authors: Savvas Paragkamian
##
## Date Created: 2025-02-01
##
## Purpose of script:
## Map generation functions for the Article 17 reporting outputs.
## Produces overview maps (Natura 2000 base, Art. 17 all species,
## border-quality QC) and per-species maps (occurrences, range,
## distribution density). All functions return the written file
## path(s) for use with targets format = "file".
##
## Funded by: Natural Environment & Climate Change Agency (NECCA/ΟΦΥΠΕΚΑ)

library(sf)
library(dplyr)
library(ggplot2)
library(ggnewscale)

.datasets_colors <- c(
  "GBIF"                 = "seagreen",
  "NECCA_redlist"        = "#B31319",
  "E1X_MDPP_2014_2024"   = "#FDF79C",
  "E1X_DB"               = "#2BA09F",
  "E2X_DB"               = "maroon1",
  "E2X_DB_references"   = "green1",
  "E1X_DB_references"    = "#141D43",
  "DistrMap_2013_2018"   = "darkorchid1",
  "Action Plan 2019"     = "#F85C29"
)

.natura_colors <- c(
  "SCI"    = "#E69F00",
  "SPA"    = "#56B4E9",
  "SCISPA" = "#CC79A7"
)

.build_base_n2000_plot <- function(greece_regions_etrs89, natura2000_etrs89) {
  ggplot2::ggplot() +
    ggplot2::geom_sf(greece_regions_etrs89, mapping = ggplot2::aes()) +
    ggplot2::geom_sf(natura2000_etrs89,
                     mapping = ggplot2::aes(fill = SITETYPE),
                     alpha = 0.3, na.rm = FALSE, show.legend = TRUE) +
    ggplot2::scale_fill_manual(
      values = .natura_colors,
      guide  = ggplot2::guide_legend(
        override.aes = list(linetype = "solid", shape = NA)
      ),
      name = "Natura2000"
    ) +
    ggplot2::theme_bw()
}

save_borders_quality_map <- function(
  species_samples_presence_border,
  greece_regions,
  eea_grid_10km,
  path
) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  greece_etrs89   <- sf::st_transform(greece_regions, 3035)
  eea_10km_etrs89 <- sf::st_transform(eea_grid_10km, 3035)

  pts <- species_samples_presence_border
  points_inside    <- pts[pts$minimumDistanceFromBorders == 0, ]
  points_lt_500m   <- pts[pts$minimumDistanceFromBorders > 0 &
                            pts$minimumDistanceFromBorders <= 500, ]
  points_gt_500m   <- pts[pts$minimumDistanceFromBorders > 500, ]

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(eea_10km_etrs89,
                     mapping = ggplot2::aes(color = "EEA grid"),
                     fill = "transparent", show.legend = TRUE) +
    ggplot2::geom_sf(greece_etrs89,
                     mapping = ggplot2::aes(color = "Greece borders"),
                     fill = NA, show.legend = TRUE) +
    ggplot2::geom_sf(points_inside,
                     mapping = ggplot2::aes(color = "Points inside"),
                     shape = 1, show.legend = TRUE) +
    ggplot2::geom_sf(points_lt_500m,
                     mapping = ggplot2::aes(color = "Points away < 500m"),
                     shape = 1, show.legend = TRUE) +
    ggplot2::geom_sf(points_gt_500m,
                     mapping = ggplot2::aes(color = "Points away > 500m"),
                     shape = 1, show.legend = TRUE) +
    ggplot2::scale_color_manual(
      name   = "Feature Type",
      values = c(
        "EEA grid"            = "grey70",
        "Greece borders"      = "blue",
        "Points inside"       = "black",
        "Points away < 500m"  = "red",
        "Points away > 500m"  = "blue"
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.87, 0.75))

  ggplot2::ggsave(path, plot = p,
                  height = 20, width = 25, dpi = 300, units = "cm", device = "png")
  path
}

save_natura_base_map <- function(greece_regions, natura2000, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  g <- .build_base_n2000_plot(
    greece_regions_etrs89 = sf::st_transform(greece_regions, 3035),
    natura2000_etrs89     = sf::st_transform(natura2000, 3035)
  )

  ggplot2::ggsave(path, plot = g,
                  height = 20, width = 20, dpi = 300, units = "cm", device = "png")
  path
}

save_art17_overview_map <- function(
  species_samples_presence_final_sf,
  greece_regions,
  natura2000,
  path
) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  presence_sf <- species_samples_presence_final_sf

  g <- .build_base_n2000_plot(
    greece_regions_etrs89 = sf::st_transform(greece_regions, 3035),
    natura2000_etrs89     = sf::st_transform(natura2000, 3035)
  ) +
    ggplot2::geom_sf(presence_sf,
                     mapping = ggplot2::aes(shape = basisOfRecord, color = collectionCode),
                     size = 2, alpha = 0.8, show.legend = TRUE) +
    ggplot2::scale_color_manual(values = .datasets_colors, name = "Datasets") +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(position = "right",
                                    override.aes = list(linetype = 0, color = NA)),
      shape = ggplot2::guide_legend(position = "right"),
      color = ggplot2::guide_legend(position = "right",
                                    override.aes = list(linetype = 0, fill = NA))
    ) +
    ggplot2::theme(legend.position = "right")

  ggplot2::ggsave(path, plot = g,
                  height = 30, width = 40, dpi = 300, units = "cm", device = "png")
  path
}

save_species_occurrence_maps <- function(
  species_samples_presence_final_sf,
  greece_regions,
  natura2000,
  maps_dir
) {
  out_dir <- file.path(maps_dir, "species_maps")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  g_base <- .build_base_n2000_plot(
    greece_regions_etrs89 = sf::st_transform(greece_regions, 3035),
    natura2000_etrs89     = sf::st_transform(natura2000, 3035)
  )

  presence_sf  <- species_samples_presence_final_sf
  species_list <- sort(unique(presence_sf$species))

  paths <- vapply(species_list, function(sp) {
    sp_occ <- dplyr::filter(presence_sf, species == sp)
    col_f  <- .datasets_colors[unique(sp_occ$collectionCode)]

    p <- g_base +
      ggplot2::geom_sf(sp_occ,
                       mapping = ggplot2::aes(shape = basisOfRecord, color = collectionCode),
                       size = 2, alpha = 0.9, show.legend = TRUE) +
      ggplot2::scale_color_manual(values = col_f, name = "Datasets") +
      ggplot2::guides(
        fill  = ggplot2::guide_legend(position = "right",
                                      override.aes = list(linetype = 0, color = NA)),
        shape = ggplot2::guide_legend(position = "right"),
        color = ggplot2::guide_legend(position = "right",
                                      override.aes = list(linetype = 0, fill = NA))
      ) +
      ggplot2::ggtitle(sp) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title            = ggplot2::element_blank(),
        axis.text             = ggplot2::element_text(colour = "black"),
        legend.title          = ggplot2::element_text(size = 8),
        legend.position       = "right",
        legend.box.background = ggplot2::element_blank()
      )

    out_path <- file.path(out_dir, paste0("map_", sp, "_occurrences.png"))
    ggplot2::ggsave(out_path, plot = p,
                    height = 28, width = 35, dpi = 300, units = "cm", device = "png")
    out_path
  }, character(1))

  unname(paths)
}

save_species_range_maps <- function(
  species_range,
  species_samples_presence_final_sf,
  greece_regions,
  natura2000,
  maps_dir
) {
  out_dir <- file.path(maps_dir, "species_maps")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  g_base <- .build_base_n2000_plot(
    greece_regions_etrs89 = sf::st_transform(greece_regions, 3035),
    natura2000_etrs89     = sf::st_transform(natura2000, 3035)
  )

  presence_sf        <- species_samples_presence_final_sf
  colors_cell_origin <- c("range" = "mediumvioletred", "distribution" = "limegreen")
  species_list       <- sort(unique(species_range$species))

  paths <- vapply(species_list, function(sp) {
    sp_occ   <- dplyr::filter(presence_sf, species == sp)
    sp_range <- dplyr::filter(species_range, species == sp)
    col_f       <- .datasets_colors[unique(sp_occ$collectionCode)]
    col_range_f <- colors_cell_origin[unique(sp_range$cell_origin)]

    p <- g_base +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_sf(sp_range,
                       mapping = ggplot2::aes(fill = cell_origin),
                       alpha = 0.8, colour = "transparent", na.rm = FALSE) +
      ggplot2::scale_fill_manual(
        values = col_range_f,
        guide  = ggplot2::guide_legend(
          override.aes = list(linetype = "solid", shape = NA)
        ),
        name = "Cell origin"
      ) +
      ggplot2::geom_sf(sp_occ,
                       mapping = ggplot2::aes(color = collectionCode),
                       size = 1, alpha = 0.9, show.legend = TRUE) +
      ggplot2::scale_color_manual(values = col_f, name = "Datasets") +
      ggplot2::guides(
        fill = ggplot2::guide_legend(position = "right",
                                     override.aes = list(linetype = 0, color = NA))
      ) +
      ggplot2::ggtitle(sp) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title            = ggplot2::element_blank(),
        axis.text             = ggplot2::element_text(colour = "black"),
        legend.title          = ggplot2::element_text(size = 8),
        legend.position       = "right",
        legend.box.background = ggplot2::element_blank()
      )

    out_path <- file.path(out_dir, paste0("map_", sp, "_range.png"))
    ggplot2::ggsave(out_path, plot = p,
                    height = 28, width = 35, dpi = 300, units = "cm", device = "png")
    out_path
  }, character(1))

  unname(paths)
}

save_species_range_shp <- function(species_range, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(species_range, path, delete_layer = TRUE, quiet = TRUE)
  path
}

save_species_distribution_maps <- function(
  species_samples_presence_final_sf,
  eea_grid_10km,
  greece_regions,
  natura2000,
  maps_dir
) {
  out_dir <- file.path(maps_dir, "species_maps")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  eea_10km_etrs89 <- sf::st_transform(eea_grid_10km, 3035)
  g_base <- .build_base_n2000_plot(
    greece_regions_etrs89 = sf::st_transform(greece_regions, 3035),
    natura2000_etrs89     = sf::st_transform(natura2000, 3035)
  )

  presence_sf  <- species_samples_presence_final_sf
  species_list <- sort(unique(presence_sf$species))

  paths <- vapply(species_list, function(sp) {
    sp_occ <- dplyr::filter(presence_sf, species == sp)
    col_f  <- .datasets_colors[unique(sp_occ$collectionCode)]

    # Grid cells coloured by number of distinct sampled locations
    locations_grid <- eea_10km_etrs89 |>
      dplyr::inner_join(
        sp_occ |>
          sf::st_drop_geometry() |>
          dplyr::distinct(CELLCODE_eea_10km, decimalLatitude, decimalLongitude) |>
          dplyr::group_by(CELLCODE_eea_10km) |>
          dplyr::summarise(n_samples = dplyr::n(), .groups = "keep") |>
          dplyr::rename(CELLCODE = CELLCODE_eea_10km),
        by = "CELLCODE"
      )

    p <- g_base +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_sf(locations_grid,
                       mapping = ggplot2::aes(fill = n_samples),
                       alpha = 0.5, colour = "transparent", na.rm = FALSE) +
      ggplot2::scale_fill_viridis_c(option = "viridis", direction = 1) +
      ggplot2::geom_sf(sp_occ,
                       mapping = ggplot2::aes(color = collectionCode),
                       size = 0.5, alpha = 0.5) +
      ggplot2::scale_color_manual(values = col_f, name = "Datasets") +
      ggplot2::guides(
        color = ggplot2::guide_legend(position = "right",
                                      override.aes = list(linetype = 0, fill = NA)),
        fill  = ggplot2::guide_colourbar(
          position = "right", alpha = 1, ticks = FALSE, label = TRUE,
          title = "n_samples", title.vjust = 0.8, order = 1
        )
      ) +
      ggplot2::ggtitle(sp) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title            = ggplot2::element_blank(),
        axis.text             = ggplot2::element_text(colour = "black"),
        legend.title          = ggplot2::element_text(size = 8),
        legend.position       = "right",
        legend.box.background = ggplot2::element_blank()
      )

    out_path <- file.path(out_dir, paste0("map_", sp, "_distribution.png"))
    ggplot2::ggsave(out_path, plot = p,
                    height = 28, width = 35, dpi = 300, units = "cm", device = "png")
    out_path
  }, character(1))

  unname(paths)
}
