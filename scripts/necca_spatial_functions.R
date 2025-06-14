library(sf)
library(dplyr)
library(magrittr)
library(terra)

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
    grids <- distribution # species distribution

    # 1. Calculate centroid distances
    centroids <- st_centroid(grids)
    dist_matrix <- st_distance(centroids)
    dist_num <- drop_units(dist_matrix)
    
    # 2. Name rows/columns of matrix by CELLCODE
    cellcode <- grids[[cellcode_col]]
    rownames(dist_num) <- colnames(dist_num) <- cellcode
    
    # 3. Create adjacency matrix: keep only upper triangle and filter by gap
    upper_only <- dist_num
    upper_only[!upper.tri(upper_only)] <- 0
    adj_matrix <- upper_only
    adj_matrix[adj_matrix <= 0 | adj_matrix > gap_distance_m] <- 0
    diag(adj_matrix) <- 0
    
    # 4. Convert adjacency matrix to edgelist
    edge_list <- which(adj_matrix != 0, arr.ind = TRUE)
    edge_df <- data.frame(
        from = rownames(adj_matrix)[edge_list[, 1]],
        to   = colnames(adj_matrix)[edge_list[, 2]],
        weight = adj_matrix[edge_list]
    )
    
    # 5. Join centroid geometries for from/to
    centroids_df <- centroids %>% st_drop_geometry() %>% mutate(geometry = st_geometry(centroids))
    from_geom <- centroids[match(edge_df$from, cellcode), ]
    to_geom   <- centroids[match(edge_df$to, cellcode), ]
    
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
    
    # 7. Find cells in full grid intersecting these lines
    intersecting_cells <- st_intersects(full_grid, lines_sf, sparse = FALSE)
    filled_cells <- full_grid[apply(intersecting_cells, 1, any), ] %>%
        mutate(cell_origin = "range")
    
    # 8. Return union of filled cells and original grids
    expanded_grid <- rbind(filled_cells, grids) %>% distinct()
    
    return(expanded_grid)
}

