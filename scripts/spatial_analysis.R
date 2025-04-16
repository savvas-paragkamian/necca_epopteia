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

########################## Load Data ##########################
###
### Species occurrences enriched ######
species_occurrences_art17_invertebrates <- read_delim("../results/species_occurrences_art17_invertebrates.tsv",delim="\t")

species_occurrences_art17_sf <- species_occurrences_art17_invertebrates |> 
    filter(!is.na(decimalLatitude)) |>
    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),
             remove=F,
             crs="WGS84")

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


### Ecosystem Types of Europe 3.1 Terrestrial
ecosystem_types_gr <- rast("/Users/talos/Documents/programming_projects/necca_epopteia/spatial_data/ecosystem_types_gr/crop_eea_r_3035_100_m_etm-terrestrial-c_2012_v3-1_r00.tif")

### Vegetation_map_Greece
## the encoding ENCODING=WINDOWS-1253 helped to see the greek characters
vegetation_map <- sf::st_read("../spatial_data/Vegetation_map_Greece/D_xabxg_VPG_60-98_GEO.shp",
                              options = "ENCODING=WINDOWS-1253")

vegetation_map_wgs <- st_transform(vegetation_map,4326)

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

EUNIS_Habitats_wgs <- st_transform(EUNIS_Habitats,4326)

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

########################### Extract data to species occurrences ############

### rasters stack

rasters_list <- list(eu_dem_gr,eu_dem_slope,ecosystem_types_gr) # hilda_rast_list, add later ,world_clim_list

points_sf <- species_occurrences_art17_sf

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


write_delim(results_ext_3,"../results/species_occurrences_art17_spatial.tsv",delim="\t")


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



######################
# Define bounding box as sf object for visual checks if needed
bbox <- st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = 4326)
bbox_sf <- st_as_sfc(bbox)

# Generate 4 random points inside the bounding box
set.seed(123)
points_inside <- st_as_sf(data.frame(
  id = 1:4,
  x = runif(4, xmin, xmax),
  y = runif(4, ymin, ymax)
), coords = c("x", "y"), crs = 4326)

# Generate 4 random points outside the bounding box
points_outside <- st_as_sf(data.frame(
  id = 5:8,
  x = c(runif(2, xmin - 10, xmin - 1), runif(2, xmax + 1, xmax + 10)),
  y = c(runif(2, ymin - 10, ymin - 1), runif(2, ymax + 1, ymax + 10))
), coords = c("x", "y"), crs = 4326)

# Combine the points into one sf object for easy handling
all_points <- rbind(points_inside, points_outside)



