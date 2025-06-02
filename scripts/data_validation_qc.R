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
library(ggnewscale)

###
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")

greece_regions_ETRS89 <- st_transform(greece_regions,3035)

N2000_v32 <- sf::st_read("../spatial_data/N2000_spatial_GR_2021_12_09_v32/N2000_spatial_GR_2021_12_09_v32.shp")

N2000_v32_ETRS89 <- st_transform(N2000_v32, 3035)

########################### Validation - QC ###########################
######### Anadohos ExtendedTemplate_Reporting dataSpecies_total ##################
#### 2025-06-02
#### load files
ranges <- sf::st_read("../anadoxos_deliverables/Maps/Range_FINAL_LANDgadm.shp")
population <- sf::st_read("../anadoxos_deliverables/Maps/Population1X1_FINAL_LANDgadm.shp")
distribution <- sf::st_read("../anadoxos_deliverables/Maps/Distribution_FINAL_LANDgadm.shp")

#### summary
#### distribution
distribution_summary <- distribution |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    st_drop_geometry() |>
    group_by(species) |>
    summarise(total_area=round(drop_units(sum(area)),digits=2))

write_delim(distribution_summary, "../anadoxos_deliverables/validation_distribution.tsv",delim="\t")

#### range summary
range_summary <- ranges |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    st_drop_geometry() |>
    group_by(Species) |>
    summarise(total_area=round(drop_units(sum(area)),digits=2))

write_delim(range_summary, "../anadoxos_deliverables/validation_range.tsv",delim="\t")
#### check distinct cells
range_dist_summary <- ranges |>
    distinct(Species,CELLCODE,geometry) |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    st_drop_geometry() |>
    group_by(Species) |>
    summarise(n_cells=n(),
              total_area=round(drop_units(sum(area)),digits=2))

#sf_use_s2(TRUE)
#### population, calculate 

population_summary <- population |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    mutate(area_tr=units::set_units(st_area(st_transform(geometry,4326)),"km^2")) |>
    st_drop_geometry() |>
    group_by(species) |>
    summarise(
              total_area=round(drop_units(sum(area)),digits=2),
              total_area_tr=round(drop_units(sum(area_tr)),digits=2))

#### population check distinct population 
population_distinct <- population |>
    distinct(species,CELLCODE_1,geometry)

#### there are duplicate cells 

population_distinct_area <- population |>
    distinct(species,CELLCODE_1,geometry) |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    mutate(area_tr=units::set_units(st_area(st_transform(geometry,4326)),"km^2")) |>
    st_drop_geometry() |>
    group_by(species) |>
    summarise(
              n_cells=n(),
              total_area=round(drop_units(sum(area)),digits=2),
              total_area_tr=round(drop_units(sum(area_tr)),digits=2))

write_delim(population_distinct_area, "../anadoxos_deliverables/validation_population.tsv",delim="\t")

#### population over N2K

population_n2k_area <- population |>
    distinct(species,CELLCODE_1,geometry)

# find the cells that overlap with N2K
intersection <- st_intersects(population_distinct, N2000_v32_ETRS89)
# find the elements with intersection
intersection_index <- lengths(intersection) > 0 
# find the elements NO with intersection
no_intersection_index <- lengths(intersection) == 0 
# keep only the intersected elements
intersection_shp <- population_distinct[intersection_index, ] 
# remove the intersected elements
no_intersection_shp <- population_distinct[no_intersection_index, ] 

# do it at once with mutate
population_distinct_n2k <- population_distinct |> 
    mutate(n2k = lengths(intersection) > 0)

nrow(intersection_shp) + nrow(no_intersection_shp)==nrow(population_distinct)


population_distinct_n2k_area <- population_distinct_n2k |>
    filter(n2k==TRUE) |>
    distinct(species,CELLCODE_1,geometry) |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    mutate(area_tr=units::set_units(st_area(st_transform(geometry,4326)),"km^2")) |>
    st_drop_geometry() |>
    group_by(species) |>
    summarise(
              n_cells=n(),
              total_area=round(drop_units(sum(area)),digits=2),
              total_area_tr=round(drop_units(sum(area_tr)),digits=2))


### validate with map


natura_colors <- c(
                   "SCI"="#E69F00",
                   "SPA"="#56B4E9",
                   "SCISPA"="#CC79A7"
)


## natura2000
g_base_n2000 <- ggplot()+
    geom_sf(greece_regions_ETRS89, mapping=aes()) +
    geom_sf(N2000_v32_ETRS89, mapping=aes(fill=SITETYPE),
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
        guides(
               fill=guide_legend(position = "inside",override.aes = list(alpha=1, linetype = 0,color=NA))
               )+
        new_scale_fill()+
        geom_sf(population_distinct_n2k, mapping=aes(fill=n2k),
                alpha=0.8,
                colour="transparent",
                na.rm = F) +
        scale_fill_manual(values=c("red","blue"),
                            guide = "fill")+
        guides(
               fill=guide_legend(position = "inside",override.aes = list(alpha=1, linetype = 0,color=NA)))+

    theme_bw() +
    theme(legend.position.inside = c(0.87, 0.75))

ggsave("../figures/map_validation_population_natura.png", 
           plot=g_base_n2000, 
           height = 20, 
           width = 20,
           dpi = 300, 
           units="cm",
           device="png")

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
