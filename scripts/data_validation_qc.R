#!/usr/bin/Rscript

## Script name: data_validation_qc.R
##
## Purpose of script:
##
## Author: Savvas Paragkamian
##
## Date Created: 2025-02-14

library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(units)

###
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")
########################### Validation - QC proigoumeni epopteia ##################

deigmata_data <- read_xlsx("../data/Invertebrates_v1.5_ALL.xlsx",
                           sheet="Δείγματα Ασπόνδυλων",
                           col_names=T
                           ) |> slice(-1)

deigmata_data_qc <- deigmata_data |> 
    mutate(latitude=as.numeric(`Γεωγραφικό Πλάτος (WGS84) Αρχη`),
           longitude=as.numeric(`Γεωγραφικό Μήκος (WGS84) Αρχή`)) |> 
    filter(!is.na(longitude))

eidi_data <- read_xlsx("../data/Invertebrates_v1.5_ALL.xlsx",
                       sheet="Είδη",
                       col_names=T
                           ) |> slice(-1)



## missing coordinates
deigmata_data_no_coords <- deigmata_data |> 
    filter(is.na('Γεωγραφικό Μήκος (WGS84) Αρχή'))


deigmata_data_sf <- st_as_sf(deigmata_data_qc,
                                   coords=c("longitude","latitude"),
                                   remove=F,
                                   crs="WGS84")


epopteia_previous_species_gr_map <- ggplot() +
    geom_sf(greece_regions, mapping=aes()) +
    geom_point(deigmata_data_qc,
            mapping=aes(x=longitude,
                        y=latitude,
                        color=`Χρηματοδοτικό Μέσο`),
            size=1.8,
            alpha=0.8,
            show.legend=T) +
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "bottom",
          legend.box.background = element_blank())

ggsave("../figures/epopteia_previous_species_gr_map.png", 
       plot=epopteia_previous_species_gr_map, 
       height = 40, 
       width = 40,
       dpi = 300, 
       units="cm",
       device="png")





###################### dadia epopteia 2015 #################################
###
#dadia <- sf::st_read("../spatial_data/Β2.Aspondyla/DeigmaAspondyla.shp")
#dadia_species <- sf::st_read("../spatial_data/Β2.Aspondyla/DeigmaAspondylaXSpecies.shp")
#
## transform 
#
#dadia_wgs <- dadia |>
#    st_transform(wgs84) 
#
#dadia_wgs_c <-  cbind(dadia_wgs, st_coordinates(dadia_wgs))
#
#dadia_d <- dadia_wgs_c |> st_drop_geometry() |> as.tibble()
### 
#dadia_s_wgs <- dadia_species |>
#    st_transform(wgs84) 
#
#dadia_s_wgs_c <-  cbind(dadia_s_wgs, st_coordinates(dadia_s_wgs))
#
#dadia_s_d <- dadia_s_wgs_c |> st_drop_geometry() |> as.tibble()
#
#
#write_delim(dadia_d, "../results/dadia_d.tsv","\t")
#write_delim(dadia_s_d, "../results/dadia_s_d.tsv","\t")
#
#st_write(dadia_wgs,"../results/dadia_wgs/dadia_wgs.shp")
#
