#!/usr/bin/Rscript

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
### greece 

greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")

#EEA 1km
eea_1km <- sf::st_read("../spatial_data/eea_reference_grid/gr_1km.shp")
eea_10km <- sf::st_read("../spatial_data/eea_reference_grid/gr_10km.shp")

### World Clim , Bioclim Variables
###
world_clim_directory <- "../spatial_data/world_clim_greece/"
world_clim_files <- list.files(world_clim_directory, pattern = "\\.tif$", full.names = TRUE)

world_clim_list <- lapply(world_clim_files, rast)

### Natura2000

N2000_v32 <- sf::st_read("../spatial_data/N2000_spatial_GR_2021_12_09_v32/N2000_spatial_GR_2021_12_09_v32.shp")

#01_Χαρτογράφηση χερσαίων Τ.Ο._EKXA





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


#### Hilda analysis
hilda_path <- "../spatial_data/hilda_greece/"
hilda_id_names <- read_delim(paste0(hilda_path, "hilda_transitions_names.tsv", sep=""), delim="\t")
hilda_all <- list.files(hilda_path)
hilda_files <- hilda_all[grepl("*.tif", hilda_all)] 


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
    
    ggsave(paste0("../figures/hilda_greece/crete_",raster_name,".png",sep=""), 
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



