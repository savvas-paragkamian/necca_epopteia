## Script name: helper_functions.R
##
## Authors: Savvas Paragkamian, Christina Kassara
##
## Date Created: 2024-04-30
##
## Purpose of script:
## Utility functions shared across the pipeline, organised in three groups:
##   1. Spatial helpers — raster cropping, polygon/raster attribute extraction,
##      range expansion with gap-distance buffering.
##   2. GBIF download helpers — resolve species names to taxon keys, submit
##      and retrieve asynchronous GBIF occurrence downloads.
##   3. Pipeline-specific raster preparation — crop large source rasters
##      (EU DEM, etc.) to Greece before first pipeline run.
##   4. Taxonomic name verification — verify species names against Catalogue
##      of Life, WoRMS, GBIF and EOL; output requires manual review before
##      feeding the pipeline as species_taxonomy_curated.tsv.
##
## Funded by: Natural Environment & Climate Change Agency (NECCA/ΟΦΥΠΕΚΑ)

library(sf)
library(dplyr)
library(magrittr)
library(terra)

## function to crop tif rasters based on bbox of other shapefile
## and make it wgs84 <- "EPSG:4326" 

#raster_path <- eurodem_file
#bbox_sf <- bbox_polygon
#output_path <- eurodem_gr_d
crop_raster_by_shp <- function(raster_path, shp, output_path) {
    # Read the raster
    raster_tmp <- rast(raster_path)
    bbox_sf <- st_as_sf(st_as_sfc(st_bbox(shp)))
    # Check CRS match
    raster_crs <- crs(raster_tmp, proj=TRUE)
    bbox_crs <- st_crs(bbox_sf)$wkt

    if (raster_crs != bbox_crs) {
        message("Reprojecting bbox to match raster CRS...")
        bbox_sf <- st_transform(bbox_sf, crs = st_crs(raster_crs))
    }
    # Convert bbox to SpatVector for terra cropping
    #bbox_vect <- vect(bbox_sf)
    #Crop raster
    cropped_raster <- crop(raster_tmp, bbox_sf)
    # project to WGS84
    wgs84 <- "EPSG:4326"
    if (crs(cropped_raster) != wgs84) {
        cropped_raster <- project(cropped_raster, wgs84)
    }
    mask_raster <- mask(cropped_raster, shp)

    output_path <- paste0(output_path,"crop_",basename(raster_path),sep="")
  
    # Save the cropped raster
    terra::writeRaster(mask_raster, output_path, overwrite = TRUE)
  
    # Clean up
    rm(raster_tmp, cropped_raster, mask_raster)
    gc()
  
    message("Saved cropped raster to: ", output_path)
}
#### extract info from polygons to add as metadata to points
extract_polygon_info_multi <- function(points_sf, polygons_list, suffixes = NULL) {
    if (!inherits(points_sf, "sf")) stop("points_sf must be an sf object.")
    if (!is.list(polygons_list)) stop("polygons_list must be a list of sf objects.")
    
    result <- points_sf
    
    for (i in seq_along(polygons_list)) {
        poly_sf <- polygons_list[[i]]
        if (!inherits(poly_sf, "sf")) stop(paste0("Item ", i, " in polygons_list is not an sf object."))
        
        # Perform spatial join
        print(polygons_list[i])
        joined <- st_join(result, poly_sf, join = st_within, left = TRUE, largest = FALSE)
        
        # Detect new columns added by join (excluding geometry)
        new_cols <- setdiff(names(joined), names(result))
        new_cols <- new_cols[new_cols != attr(joined, "sf_column")]  # exclude geometry column
        
        # If no new columns, skip this layer
        if (length(new_cols) == 0) {
          warning(paste("No new columns added from polygon layer", i, "- skipping."))
          next
        }
  
        # Add suffixes to new columns
        suffix <- if (!is.null(suffixes) && length(suffixes) >= i) suffixes[i] else paste0("poly", i)
        names(joined)[names(joined) %in% new_cols] <- paste0(new_cols, "_", suffix)
        
        # Drop geometry to avoid duplication
        joined_no_geom <- joined %>% st_drop_geometry()
        
        # Append only the new columns
        result <- bind_cols(result, joined_no_geom[, paste0(new_cols, "_", suffix), drop = FALSE])
    }
    
    return(result)
}
#### extract info from rasters to add to points as metadata
extract_from_named_rasters <- function(raster_list, points_sf) {
    if (!all(sapply(raster_list, inherits, "SpatRaster"))) {
        stop("All items in raster_list must be terra SpatRaster objects.")
    }
  
    # Get names of each raster in the list
    raster_names <- names(raster_list)
  
    # Convert sf to terra::vect
    vect_points <- terra::vect(points_sf)
  
    # Loop through rasters and extract
    extracted_list <- lapply(seq_along(raster_list), function(i) {
        r <- raster_list[[i]]
        vals <- terra::extract(r, vect_points)[, -1, drop = FALSE]
  
        # Generate nice column names using object name and band names
        layer_names <- names(r)
        if (is.null(layer_names) || any(layer_names == "")) {
          layer_names <- paste0("band", seq_len(ncol(vals)))
        }
  
        prefix <- raster_names[i]
        colnames(vals) <- paste0(prefix, "_", layer_names)
  
        return(vals)
    })
  
    # Combine all extracted values
    all_extracted <- do.call(cbind, extracted_list)
  
    # Add to sf object
    result_sf <- cbind(points_sf, all_extracted)
  
    return(result_sf)
}


### range using the centroid

expand_range_with_gap_distance <- function(distribution, full_grid, gap_distance_m = 40000, cellcode_col = "CELLCODE") {

    # 0. initiate
    grids <- distribution |>
        mutate(cell_origin="distribution")# species distribution

    # 1. Calculate centroid distances
    centroids <- st_centroid(grids)
    dist_matrix <- st_distance(centroids, which="Euclidean")
    dist_num <- drop_units(dist_matrix)
    #print("centroids")
    
    # 2. Name rows/columns of matrix by CELLCODE
    cellcode <- grids[[cellcode_col]]
    rownames(dist_num) <- colnames(dist_num) <- cellcode
    
    # 3. Create adjacency matrix: keep only upper triangle and filter by gap
    upper_only <- dist_num
    upper_only[!upper.tri(upper_only)] <- 0
    adj_matrix <- upper_only
    adj_matrix[adj_matrix <= 0 | adj_matrix > gap_distance_m] <- 0
    diag(adj_matrix) <- 0
    
    #print("matrix")
    # 4. Convert adjacency matrix to edgelist
    edge_list <- which(adj_matrix != 0, arr.ind = TRUE)
    edge_df <- data.frame(
        from = rownames(adj_matrix)[edge_list[, 1]],
        to   = colnames(adj_matrix)[edge_list[, 2]],
        weight = adj_matrix[edge_list]
    )
    
    #print("edgelist")
    # 5. Join centroid geometries for from/to
    centroids_df <- centroids %>% st_drop_geometry() %>% mutate(geometry = st_geometry(centroids))
    from_geom <- centroids[match(edge_df$from, cellcode), ]
    to_geom   <- centroids[match(edge_df$to, cellcode), ]
    
    #print("join centroid")
    # 6. Create lines between centroids
    line_list <- mapply(
      function(p1, p2) st_linestring(matrix(c(st_coordinates(p1), st_coordinates(p2)), ncol = 2, byrow = TRUE)),
      st_geometry(from_geom),
      st_geometry(to_geom),
      SIMPLIFY = FALSE
    )
    
    lines_sf <- st_sf(
        from = edge_df$from,
        to   = edge_df$to,
        geometry = st_sfc(line_list, crs = st_crs(grids))
    )
    
    #print("connect with line")
    # 7. Find cells in full grid intersecting these lines
    intersecting_cells <- st_crosses(full_grid, lines_sf, sparse = FALSE)
    filled_cells <- full_grid[apply(intersecting_cells, 1, any), ] %>%
        mutate(cell_origin = "range")
    #print(colnames(filled_cells))
    #print(colnames(grids))
    
    #print("find grids intersect with lines")
    # 8. Return union of filled cells and original grids
    expanded_grid <- rbind(filled_cells, grids) %>% distinct()
    
    print("range calculated")
    return(expanded_grid)
}

# =============================================================================
# GBIF download helpers
# These functions are interactive / one-time preparation tools.
# GBIF credentials must be set in ~/.Renviron:
#   GBIF_USER, GBIF_PWD, GBIF_EMAIL
# =============================================================================

# Resolve a character vector of species names to GBIF backbone taxon keys.
# Returns a numeric vector of unique usageKey values.
get_gbif_taxon_keys <- function(species_names) {
  result <- rgbif::name_backbone_checklist(species_names, verbose = TRUE)
  keys <- result |>
    dplyr::filter(is_alternative == "FALSE", !is.na(usageKey)) |>
    dplyr::pull(usageKey) |>
    unique()
  message("Resolved ", length(keys), " unique GBIF taxon keys from ",
          length(species_names), " species names.")
  keys
}

# Submit an asynchronous GBIF occurrence download request.
# Returns the download key (character) needed for import_gbif_download().
# Predicates: country filter, coordinates required, no geospatial issues.
request_gbif_download <- function(taxon_keys, country = "GR") {
  dl <- rgbif::occ_download(
    rgbif::pred_in("taxonKey", as.integer(taxon_keys)),
    rgbif::pred("country", country),
    rgbif::pred("hasCoordinate", TRUE),
    rgbif::pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV"
  )
  key <- as.character(dl)
  message("Download requested. Key: ", key,
          "\nMonitor at: https://www.gbif.org/occurrence/download/", key)
  key
}

# Wait for a GBIF download to complete, import it, and save as TSV.
# download_key : the key returned by request_gbif_download()
# output_path  : destination .tsv file (directory created if absent)
# Returns output_path invisibly.
import_gbif_download <- function(download_key, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  message("Waiting for GBIF download ", download_key, " …")
  rgbif::occ_download_wait(download_key)
  occ <- rgbif::occ_download_get(download_key) |>
    rgbif::occ_download_import()
  readr::write_delim(occ, output_path, delim = "\t")
  message("Saved ", nrow(occ), " GBIF records to: ", output_path)
  invisible(output_path)
}

# Convenience wrapper: resolve names → request download → wait → save.
# Runs the full workflow in one call. Not suitable for large downloads
# (no progress feedback between request and completion).
download_gbif_occurrences <- function(species_names, output_path, country = "GR") {
  keys <- get_gbif_taxon_keys(species_names)
  key  <- request_gbif_download(keys, country = country)
  import_gbif_download(key, output_path)
}

# =============================================================================
# Raster cropping helpers
# =============================================================================

# Crop a large raster file to the bounding box of a reference sf object,
# optionally reprojecting to a target CRS before saving.
#
# raster_path : path to the source raster (any format readable by terra)
# extent_sf   : an sf object whose bounding box defines the crop extent
# output_path : destination file path (.tif); directory created if absent
# target_crs  : EPSG code (integer) to reproject the output, or NULL to keep
#               the raster's original CRS (default: NULL)
# mask        : if TRUE, mask cells outside the bounding box (default: FALSE)
#
# Returns output_path invisibly.
crop_raster_to_extent <- function(raster_path, extent_sf,
                                  output_path, target_crs = NULL,
                                  mask = FALSE) {
  r <- terra::rast(raster_path)

  bbox_sf <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(extent_sf)))
  bbox_reproj <- sf::st_transform(bbox_sf, sf::st_crs(terra::crs(r)))
  bbox_vect <- terra::vect(bbox_reproj)

  cropped <- terra::crop(r, bbox_vect)

  if (mask) {
    cropped <- terra::mask(cropped, bbox_vect)
  }

  if (!is.null(target_crs)) {
    cropped <- terra::project(cropped, paste0("EPSG:", target_crs))
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(cropped, output_path, overwrite = TRUE)
  message("Saved cropped raster to: ", output_path)
  invisible(output_path)
}

# Crop the EU DEM (stored in EPSG:3035) to Greece and save.
# Thin wrapper around crop_raster_to_extent() with pipeline-specific defaults:
# the DEM stays in LAEA (3035) — no reprojection.
#
# eu_dem_path     : path to the full EU DEM mosaic (.tif, ~20 GB)
# greece_regions  : sf object of Greece administrative boundaries
# output_path     : destination .tif in data/spatial/EU_DEM_mosaic_5deg_gr/
crop_eu_dem_to_greece <- function(eu_dem_path, greece_regions, output_path) {
  greece_laea <- sf::st_transform(greece_regions, 3035)
  crop_raster_to_extent(
    raster_path = eu_dem_path,
    extent_sf   = greece_laea,
    output_path = output_path,
    target_crs  = NULL,
    mask        = FALSE
  )
}

# =============================================================================
# Taxonomic name verification
# Interactive / one-time preparation tool.
# Output requires MANUAL REVIEW before it can feed the pipeline:
#   results → human curation → data/raw/species_taxonomy_curated.tsv
# Data sources: 1 = Catalogue of Life, 9 = WoRMS, 11 = GBIF, 12 = EOL
# =============================================================================

# Verify a character vector of species names against four authoritative
# taxonomic databases and save the results as a TSV for manual curation.
#
# species_names : character vector of scientific names to verify
# output_path   : destination .tsv file (directory created if absent)
#
# Returns the verification data frame invisibly.
# NOTE: the saved TSV is an intermediate review file, not a pipeline input.
# After inspecting the matches, manually produce species_taxonomy_curated.tsv.
verify_species_taxonomy <- function(species_names,
                                    output_path = "results/gnr_species_verifier.tsv") {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  result <- taxize::gna_verifier(
    species_names,
    all_matches  = TRUE,
    data_sources = c(1, 9, 11, 12)
  )

  readr::write_delim(result, output_path, delim = "\t")
  message(
    "Verified ", length(species_names), " names against ",
    "Catalogue of Life, WoRMS, GBIF, and EOL.\n",
    "Results saved to: ", output_path, "\n",
    "ACTION REQUIRED: review matches and update species_taxonomy_curated.tsv."
  )

  invisible(result)
}
