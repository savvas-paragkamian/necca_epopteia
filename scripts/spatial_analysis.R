#!/usr/bin/env Rscript

## Script name: spatial_analysis.R
##
## Purpose of script:
##
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-11-06

library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(units)
library(ggpubr)
library(ggnewscale)

########################## Load Data ##########################
###
### Species occurrences enriched ######
species_samples_art17 <- read_delim("../results/species_samples_art17.tsv", delim="\t")
species_occurrences_invertebrates <- read_delim("../results/species_occurrences_invertebrates.tsv",delim="\t")

points_sf <- species_samples_art17 |> 
    filter(!is.na(decimalLatitude)) |>
    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),
             remove=F,
             crs="WGS84")

## distribution data from previous report
species_dist_national_rep <- sf::st_read("../spatial_data/National report_2013_2018_shp/GR_Art17_species_distribution.shp")
species_dist_national_rep_sens <- sf::st_read("../spatial_data/National report_2013_2018_shp/GR_Art17_species_distribution_sensitive.shp")
### greece 

greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")

### EEA reference grid
eea_1km <- sf::st_read("../spatial_data/eea_reference_grid/gr_1km.shp")
eea_1km_wgs <- st_transform(eea_1km,4326)

eea_10km <- sf::st_read("../spatial_data/eea_reference_grid/gr_10km.shp")
eea_10km_wgs <- st_transform(eea_10km,4326)
### EU DEM Greece
eu_dem_gr <- rast("../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_4326_gr.tif")

### EU DEM Greece slope
eu_dem_slope <- rast("../spatial_data/EU_DEM_slope_gr/crop_eudem_slop_3035_europe.tif")

### World Clim , Bioclim Variables
###
world_clim_directory <- "../spatial_data/world_clim_greece/"
world_clim_files <- list.files(world_clim_directory, pattern = "\\.tif$", full.names = TRUE)

world_clim_list <- lapply(world_clim_files, rast)

### Natura2000

N2000_v32 <- sf::st_read("../spatial_data/N2000_spatial_GR_2021_12_09_v32/N2000_spatial_GR_2021_12_09_v32.shp")

N2000_v32_wgs <- st_transform(N2000_v32,4326)
N2000_v32_wgs_sci <- N2000_v32_wgs |>
    filter(SITETYPE=="SCI")

N2000_v32_wgs_spa <- N2000_v32_wgs |>
    filter(SITETYPE=="SPA")

N2000_v32_wgs_scispa <- N2000_v32_wgs |>
    filter(SITETYPE=="SCISPA")


### Ecosystem Types of Europe 3.1 Terrestrial
ecosystem_types_gr <- rast("/Users/talos/Documents/programming_projects/necca_epopteia/spatial_data/ecosystem_types_gr/crop_eea_r_3035_100_m_etm-terrestrial-c_2012_v3-1_r00.tif")

### Vegetation_map_Greece
## the encoding ENCODING=WINDOWS-1253 helped to see the greek characters
vegetation_map <- sf::st_read("../spatial_data/Vegetation_map_Greece/D_xabxg_VPG_60-98_GEO.shp",
                              options = "ENCODING=WINDOWS-1253")

vegetation_map_wgs <- st_transform(vegetation_map,4326) |>
    st_make_valid()

veg_data <- data.frame(
  A_VEG_TYPE = c("ΕΛΑ", "ΕΡΛ", "ΠΜΑ", "ΠΛΔ", "ΠΔΑ", "ΠΧΑ", "ΠΚΟ", "ΠΘΑ", "ΚΠΡ",
                 "ΑΡΚ", "ΟΞΥ", "ΔΡΥ", "ΚΑΣ", "ΣΗΜ", "ΣΦΕ", "ΦΙΛ", "ΦΠΛ", "ΠΑΡ",
                 "ΕΥΚ", "ΦΟΙ", "ΘΑΜ", "ΦΘΑ", "ΛΙΒ", "ΑΓΟ", "ΟΙΚ", "ΓΚΑ", "ΓΚΕ",
                 "ΛΧΡ", "ΛΙΜ"),
  A_VEG_NAME = c("Ελάτη", "Ερυθρελάτη", "Πεύκη μαύρη", "Πεύκη λευκόδερμη", "Πεύκη δασική",
                 "Πεύκη χαλέπιος", "Πεύκη κουκουναριά", "Πεύκη θαλασσία", "Κυπαρίσσι",
                 "Άρκευθος", "Οξυά", "Δρύς", "Καστανιά", "Σημύδα", "Σφένδαμος",
                 "Φιλύρα", "Φυλλοβόλα πλατύφυλλα", "Παραποτάμια βλάστηση", "Ευκάλυπτος",
                 "Φοίνικας", "Θάμνοι", "Φυλλοβόλοι θάμνοι", "Λιβάδια, αραιά ξυλ. βλάστηση",
                 "Άγονα", "Οικισμοί", "Γεωργ. καλλιέργειες", "Γεωργ. καλλιέργειες εγκατ.",
                 "Λοιπές χρήσεις", "Λίμνη"),
  stringsAsFactors = FALSE
)

### EUNIS_Habitats_2018 for the Natura2000 areas
# the same as 01_Χαρτογράφηση χερσαίων Τ.Ο._EKXA

EUNIS_Habitats <- sf::st_read("../spatial_data/EUNIS_Habitats_2018/Habitats_2018/Habitats.shp")

EUNIS_Habitats_wgs <- st_transform(EUNIS_Habitats,4326) |>
    st_make_valid()

######### HILDA Greece land use change ##########
hilda_cat <- data.frame(hilda_id = c("11","22","33","44","55","66","77"),
                        hilda_name=c("urban","cropland","pasture/rangeland",
                                     "forest", "unmanaged grass/shrubland","sparse/no vegetation", "water"),
                        hilda_hex=c("#000000","#AE6120","#98BA6A","#07A07D","#BE81A3","#999999", "#1370A1"))
hilda_cat_v <- c("urban"="#000000",
                 "cropland"="#AE6120",
                 "pasture/rangeland"="#98BA6A",
                 "forest"="#07A07D",
                 "unmanaged grass/shrubland"="#BE81A3",
                 "sparse/no vegetation"="#999999",
                 "water"="#1370A1")


#### Hilda land cover change
hilda_path <- "../spatial_data/hilda_greece/"
hilda_id_names <- read_delim(paste0(hilda_path, "hilda_transitions_names.tsv", sep=""), delim="\t")
hilda_all <- list.files(hilda_path)
hilda_files <- hilda_all[grepl("*.tif", hilda_all)] 
hilda_rast_list <- lapply(hilda_files, function(x) rast(paste0(hilda_path,x,sep="")))

######################## Extract data to species occurrences ####################

### rasters stack

rasters_list <- list(eu_dem_gr,eu_dem_slope,ecosystem_types_gr) # hilda_rast_list, add later ,world_clim_list


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


results_ext <- extract_from_named_rasters(rasters_list,points_sf)
results_ext_2 <- extract_from_named_rasters(world_clim_list,results_ext)
results_ext_3 <- extract_from_named_rasters(hilda_rast_list,results_ext_2)


### shapefiles

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

result <- extract_polygon_info_multi(
  results_ext_3,
  polygons_list <- list(eea_1km_wgs,
                      eea_10km_wgs,
                      N2000_v32_wgs_sci,
                      N2000_v32_wgs_spa,
                      N2000_v32_wgs_scispa,
                      vegetation_map_wgs,
                      EUNIS_Habitats_wgs),
  suffixes = c("eea_1km","eea_10km","N2000_v32_sci","N2000_v32_spa","N2000_v32_scispa","vegetation_map","EUNIS_Habitats")
)



species_occurrences_spatial <- result
write_delim(species_occurrences_spatial,"../results/species_occurrences_spatial.tsv",delim="\t")


################################## MAPS ########################################
################################################################################
################################################################################

species_with_data <- unique(points_sf$species)

datasets_colors <- c(
                     "GBIF"="seagreen",
                     "NECCA_redlist"="#B31319",
                     "E1X_MDPP_2014_2024"="#FDF79C",
                     "E1X_DB"="#2BA09F",
                     "E1X_DB_references"="#141D43",
                     "Invertebrates_records_Olga"="#F85C29"
                     )

# base plot with Natura2000 areas of Greece
natura_colors <- c(
                   "SCI"="#E69F00",
                   "SPA"="#56B4E9",
                   "SCISPA"="#CC79A7"
)


## natura2000
g_base_n2000 <- ggplot()+
    geom_sf(greece_regions, mapping=aes()) +
    geom_sf(N2000_v32_wgs, mapping=aes(fill=SITETYPE),
            alpha=0.3,
            #colour="transparent",
            na.rm = F,
            show.legend=T) +
    scale_fill_manual(
                      values= natura_colors,
                       guide = guide_legend(
                                            override.aes = list(
                                                                linetype="solid",
                                                                shape = NA)
                                            ),
                       name="Natura2000"
                       )+
    theme_bw()

ggsave("../figures/map_natura.png", 
           plot=g_base_n2000, 
           height = 20, 
           width = 20,
           dpi = 300, 
           units="cm",
           device="png")

### natura2000 with all points
g_art17_n2000 <- g_base_n2000 +
    geom_point(points_sf,
            mapping=aes(x=decimalLongitude,
                        y=decimalLatitude,
                        color=datasetName),
            size=1.2,
            alpha=0.8,
            show.legend=T) +
    scale_color_manual(values=datasets_colors,
                        name = "Datasets")+
    guides(
           fill=guide_legend(position = "inside",override.aes = list(linetype = 0,color=NA)),
           color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)))+
    theme(legend.position.inside = c(0.87, 0.75)
    )


ggsave("../figures/map_art17_invertebrates_natura.png", 
           plot=g_art17_n2000, 
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")


############################## Species Occurrences ##############################

### figures of each invertebrate of art17 for Greece

for (i in seq_along(species_with_data)){
    species_occurrences <- points_sf |>
        filter(species==species_with_data[i])

    dataset_colors_f <- datasets_colors[unique(species_occurrences$datasetName)]
    print(species_with_data[i])

    
    species_gr_map <- g_base_n2000 +
        geom_point(species_occurrences,
                mapping=aes(x=decimalLongitude,
                            y=decimalLatitude,
                            color=datasetName),
                size=1.8,
                alpha=0.9,
                show.legend=T) +
        coord_sf(crs="WGS84") +
        scale_color_manual(values=dataset_colors_f,
                        name = "Datasets") +
        guides(
               fill=guide_legend(position = "inside",override.aes = list(linetype = 0,color=NA)),
               color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)))+
        ggtitle(paste(species_with_data[i]))+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position.inside = c(0.87, 0.75),
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_occurrences.png", sep=""), 
           plot=species_gr_map, 
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")

}




############################## Species Population ##############################
### population is defined here as the number of 1X1km cells
### remove GBIF

for (i in seq_along(species_with_data)){
    # use the 1X1 as a proxy of population

    species_occurrences <- points_sf |>
        filter(species==species_with_data[i]) |>
        filter(datasetName!="GBIF")
    
    locations_1_grid_samples <- st_join(eea_1km_wgs, species_occurrences, left=F) |>
        distinct(geometry,CELLCODE, decimalLatitude, decimalLongitude) |>
        group_by(geometry,CELLCODE) |>
        summarise(n_samples=n(),.groups="keep")

    dataset_colors_f <- datasets_colors[unique(species_occurrences$datasetName)]

    print(species_with_data[i])

    
    species_gr_map <- g_base_n2000 +
#        geom_point(species_occurrences,
#                mapping=aes(x=decimalLongitude,
#                            y=decimalLatitude,
#                            color=datasetName),
#                size=1.8,
#                alpha=0.9,
#                show.legend=T) +
#        scale_color_manual(values=dataset_colors_f,
#                            name = "Datasets")+
        new_scale_fill()+
        geom_sf(locations_1_grid_samples, mapping=aes(fill=n_samples),
                alpha=0.8,
                colour="transparent",
                na.rm = F,
                show.legend=T) +
        scale_fill_gradient(low="gray50",
                            high="gray1",
                            guide = "colourbar")+
        guides(
               color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)),
               fill = guide_colourbar(position = "inside",
                                      alpha = 1,
                                      ticks = F,
                                      label = T,
                                      title="n_samples",
                                      title.vjust = 0.8,
                                      order = 1))+
        theme_bw()+
        ggtitle(paste(species_with_data[i]))+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position.inside = c(0.87,0.65),
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_population.png", sep=""), 
           plot=species_gr_map, 
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")

}


############################## Species Distribution ##############################

for (i in seq_along(species_with_data)){
    # use the 10X10 as a proxy of distribution

    species_occurrences <- points_sf |>
        filter(species==species_with_data[i]) |>
        filter(datasetName!="GBIF")
    
    locations_10_grid_samples <- st_join(eea_10km_wgs, species_occurrences, left=F) |>
        distinct(geometry,CELLCODE, decimalLatitude, decimalLongitude) |>
        group_by(geometry,CELLCODE) |>
        summarise(n_samples=n(),.groups="keep")

    dataset_colors_f <- datasets_colors[unique(species_occurrences$datasetName)]

    print(species_with_data[i])

    
    species_gr_map <- g_base_n2000 +
        geom_point(species_occurrences,
                mapping=aes(x=decimalLongitude,
                            y=decimalLatitude,
                            color=datasetName),
                size=1.8,
                alpha=0.9,
                show.legend=T) +
        scale_color_manual(values=dataset_colors_f,
                            name = "Datasets")+
        new_scale_fill()+
        geom_sf(locations_10_grid_samples, mapping=aes(fill=n_samples),
                alpha=0.8,
                colour="transparent",
                na.rm = F,
                show.legend=T) +
        scale_fill_gradient(low="gray50",
                            high="gray1",
                            guide = "colourbar")+
        guides(
               color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)),
               fill = guide_colourbar(position = "inside",
                                      alpha = 1,
                                      ticks = F,
                                      label = T,
                                      title="n_samples",
                                      title.vjust = 0.8,
                                      order = 1))+
        theme_bw()+
        ggtitle(paste(species_with_data[i]))+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position.inside = c(0.87,0.65),
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_distribution.png", sep=""), 
           plot=species_gr_map, 
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")

}

################################ Species Range ##############################
library(igraph)

i=1
species_occurrences <- points_sf |>
    filter(species==species_with_data[i]) |>
    filter(datasetName!="GBIF")

locations_10_grid_samples <- st_join(eea_10km_wgs, species_occurrences, left=F) |>
    distinct(geometry,CELLCODE, decimalLatitude, decimalLongitude) |>
    group_by(geometry,CELLCODE) |>
    summarise(n_samples=n(),.groups="keep")

grids <- eea_10km_wgs |>
    filter(CELLCODE %in% locations_10_grid_samples$CELLCODE) |>
    mutate(cell_origin="distribution")

centroids <- st_centroid(grids)
dist_matrix <- st_distance(centroids)
dist_num <- drop_units(dist_matrix) # remove units

# keep the cell code of the grid
cellcode <- grids$CELLCODE

rownames(dist_num) <- colnames(dist_num) <- cellcode

# STEP 2: Convert distance matrix to a graph where edges are < gap threshold
gap_distance_u <- set_units(40000, "m")   # e.g., 40 km
gap_distance <- 40000   # e.g., 40 km
upper_only <- dist_num


# Zero out lower triangle and diagonal because is symmetric
upper_only[!upper.tri(upper_only)] <- 0

adj_matrix <- upper_only  # copy
# filter only the cells with shorter gap distance
adj_matrix[adj_matrix <= 0 | adj_matrix > gap_distance] <- 0
diag(adj_matrix) <- 0  # Remove self-links

# transform the matrix to edgelsit
edge_list <- which(adj_matrix != 0, arr.ind = TRUE)

edge_df <- data.frame(from = rownames(adj_matrix)[edge_list[, 1]],
                      to   = colnames(adj_matrix)[edge_list[, 2]],
                      weight = adj_matrix[edge_list])

# make geometry
# Join centroids to edgelist (twice: once for 'from', once for 'to')
line_list <- mapply(
                    function(p1, p2) st_linestring(matrix(c(st_coordinates(p1), st_coordinates(p2)), ncol = 2, byrow = TRUE)),
                    st_geometry(edges$from_geom),
                    st_geometry(edges$to_geom),
                    SIMPLIFY = FALSE
)

# Convert list to sfc and make an sf object
lines_sf <- st_sf(
                  from = edges$from,
                  to = edges$to,
                  geometry = st_sfc(line_list, crs = st_crs(grids))
                  )

# Intersect lines with full grid to get "in-between" cells
intersecting_cells <- st_intersects(eea_10km_wgs, lines_sf, sparse = FALSE)

# 7. Keep only grid cells that intersect at least one line
filled_cells <- eea_10km_wgs[apply(intersecting_cells, 1, any), ] |>
    mutate(cell_origin="range")

# Create polygons by merging connected grids
expanded_grid <- rbind(filled_cells, grids) %>% distinct()

# Maps
species_gr_map <- ggplot()+
    geom_sf(greece_regions, mapping=aes()) +
    geom_sf(N2000_v32_wgs, mapping=aes(fill=SITETYPE),
            alpha=0.3,
            #colour="transparent",
            na.rm = F,
            show.legend=T) +
    scale_fill_manual(
                      values= natura_colors,
                      guide = guide_legend(
                                            override.aes = list(
                                                                linetype="solid",
                                                                shape = NA)
                                            ),
                       name="Natura2000"
                       )+
    geom_point(species_occurrences,
            mapping=aes(x=decimalLongitude,
                        y=decimalLatitude,
                        color=datasetName),
            size=1.8,
            alpha=0.9,
            show.legend=T) +
    scale_color_manual(values=dataset_colors_f,
                        name = "Datasets")+
    new_scale_fill()+
    geom_sf(expanded_grid, mapping=aes(fill=cell_origin),
            alpha=0.8,
            colour="transparent",
            na.rm = F,
            show.legend=T) +
    guides(
           fill=guide_legend(position = "inside",override.aes = list(linetype = 0,color=NA)),
           color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)))+
    ggtitle(paste(species_with_data[i]))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position.inside = c(0.87, 0.75),
          legend.box.background = element_blank())

ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_range.png", sep=""), 
       plot=species_gr_map, 
       height = 20, 
       width = 25,
       dpi = 300, 
       units="cm",
       device="png")


############################## Hilda analysis ##############################
for (i in 1:length(hilda_files)){
    print(i)
    filename <- hilda_files[i]
    
    raster_file <- rast(paste0(hilda_path,hilda_files[i],sep=""))
    ### terra (or raster in R) gives you very small floating point numbers instead of nice, clean integers or just zeros.
    raster_file <- round(raster_file)
    raster_name <- paste0("hilda_",gsub(".*([0-9]{4}).*", "\\1", filename),sep="")
    
    raster_df <- terra::as.data.frame(raster_file, xy=TRUE, cells=TRUE)

    raster_df <- raster_df |>
        mutate(hilda_id=round(raster_df[,4], digits = 0)) |>
        filter(hilda_id>0) |>
        mutate(hilda_id=as.character(hilda_id)) |>
        left_join(hilda_cat, by=c("hilda_id"="hilda_id"))
    
    raster_df$hilda_name <- factor(raster_df$hilda_name, levels=as.character(unique(sort(raster_df$hilda_name))))
    
    g_hilda_map <- ggplot() +
        geom_sf(greece_regions, mapping=aes()) +
        geom_raster(raster_df,
                    mapping=aes(x=x, y=y, fill=hilda_name)) +
        scale_fill_manual(values=hilda_cat_v) +
        guides(fill = guide_legend(nrow = 1)) +
        ggtitle(raster_name)+
        theme(axis.title=element_blank(),
              legend.position="bottom",
              legend.key.size = unit(4, "mm"), 
              legend.text=element_text(size=8),
              legend.title=element_blank())
    
#    ggsave(paste0("../plots/hilda_crete/crete_",raster_name,"_map.png",sep=""),
#           plot=g_hilda_map,
#           height = 10, 
#           width = 20,
#           dpi = 300, 
#           units="cm",
#           device="png")
    
    hilda_sum <- zonal(cellSize(raster_file), raster_file, "sum") |> 
        mutate(area_m2=units::set_units(area,m^2)) |>
        mutate(area=units::set_units(area/10^6, km^2)) 
    
    hilda_sum <- hilda_sum |>
        mutate(hilda_id=as.character(hilda_sum[,1])) |>
        filter(hilda_sum[,1]>0) |>
        left_join(hilda_cat) 

    land_cover_colors <- c(
                           "urban" = "#000000",
                           "cropland" = "#AE6120",
                           "pasture/rangeland" = "#98BA6A",
                           "forest" = "#07A07D",
                           "unmanaged grass/shrubland" = "#BE81A3",
                           "sparse/no vegetation" = "#999999",
                           "water" = "#1370A1"
    )
    
    hilda_sum_g <- ggplot()+
        geom_col(hilda_sum,
                 mapping= aes(y=as.numeric(area),
                              x="",
                              fill = hilda_name),
                 position = position_stack(),
                 width = 0.2) +
        scale_fill_manual(values=land_cover_colors) +
        scale_x_discrete(expand = expansion(add=c(0,0)))+
#        scale_y_continuous(breaks=seq(0,9000,1000),
#                           limits=c(0,8900),
#                           expand = c(0,0))+
        ylab("Area sq. km") +
        xlab("") +
        theme_bw()+
        theme(legend.position='none',
              axis.ticks.x=element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(), #remove major gridlines
              panel.grid.minor = element_blank()) #remove minor gridlines
    
#    ggsave(paste0("../plots/hilda_crete/crete_",raster_name,"_bar.png",sep=""),
#           plot=hilda_sum_g,
#           height = 10, 
#           width = 10,
#           dpi = 300, 
#           units="cm",
#           device="png")

    fig_hilda <- ggarrange(g_hilda_map,hilda_sum_g,
              labels = c("A", "B"),
              ncol = 2,
              nrow = 1,
              widths = c(0.85,0.15),
              font.label=list(color="black",size=15),
              common.legend = TRUE,
              legend="bottom") + bgcolor("white")
    
    ggsave(paste0("../figures/hilda_greece/greece_",raster_name,".png",sep=""), 
           plot=fig_hilda, 
           height = 10, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")
}
