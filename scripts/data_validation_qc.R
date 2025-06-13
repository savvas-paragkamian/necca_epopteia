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
source("necca_spatial_functions.R")

###
#points_final <- st_read("../results/species_samples_art17.gpkg")

wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")

greece_regions_ETRS89 <- st_transform(greece_regions,3035)

eea_1km <- sf::st_read("../spatial_data/eea_reference_grid/gr_1km.shp")
eea_1km_ETRS89 <- st_transform(eea_1km, 3035)

eea_10km <- sf::st_read("../spatial_data/eea_reference_grid/gr_10km.shp")
eea_10km_ETRS89 <- st_transform(eea_10km, 3035)

## do the st_union because greece_regions_ETRS89 is a multipolygon!!
eea_10km_ETRS89_gr <- st_intersection(eea_10km_ETRS89, st_union(greece_regions_ETRS89)) |>
    distinct(CELLCODE,geometry) |> 
    mutate(area=units::set_units(st_area(geometry),"km^2"))

N2000_v32 <- sf::st_read("../spatial_data/N2000_spatial_GR_2021_12_09_v32/N2000_spatial_GR_2021_12_09_v32.shp")

N2000_v32_ETRS89 <- st_transform(N2000_v32, 3035)

########################### Validation - QC ###########################
######### Anadohos ExtendedTemplate_Reporting dataSpecies_total ##################
#### 2025-06-02
#### load files

points <- sf::st_read("../anadoxos_deliverables/Distribution&RangeMaps_20250611/VerifiedOccurrenceDB_PlusOrphans_LAEA.shp")

ranges <- sf::st_read("../anadoxos_deliverables/Distribution&RangeMaps_20250611/Species_range_invertebrates.shp") |>
        mutate(file="range")
distribution <- sf::st_read("../anadoxos_deliverables/Distribution&RangeMaps_20250611/Species_distribution_invertebrates.shp") |>
    mutate(file="distribution")

population <- sf::st_read("../anadoxos_deliverables/Distribution&RangeMaps_20250611/National_gr1km_Pops.shp")

population_n2k <- sf::st_read("../anadoxos_deliverables/Distribution&RangeMaps_20250611/N2K_gr1km_Pops.shp")

######### points quality control

points <- points |>
    mutate(occurrence_id=row_number())
points_no_points <- points[which(is.na(points$datasetNam)),] 

length(points_no_points)


## points unique
points_u <- points |> 
#    st_drop_geometry() |>
    count(Species, decimalLat,decimalLon, geometry)

## these are of Parnassius apollo
unique(points_no_points$Species)


## my population
##
my_population <- st_join(points_u,
                             eea_1km_ETRS89,
                             join = st_intersects,
                             left = TRUE,
                             largest = FALSE)


my_population_unique <- my_population |> 
    distinct(Species,CELLCODE)

my_population_summary <- my_population_unique |>
    group_by(Species) |>
    summarize(my_n_eea_1km=n())


## my distribution
my_distribution <- st_join(points_u,
                             eea_10km_ETRS89,
                             join = st_intersects,
                             left = TRUE,
                             largest = FALSE)

my_distribution_unique <- my_distribution |>
    distinct(Species,CELLCODE) |>
    left_join(eea_10km_ETRS89_gr, by=c("CELLCODE"="CELLCODE"))

my_distribution_whole_grid <- my_distribution |>
    distinct(Species,CELLCODE) |>
    left_join(eea_10km_ETRS89, by=c("CELLCODE"="CELLCODE")) |> 
    sf::st_as_sf()

my_distribution_summary <- my_distribution_unique |>
    group_by(Species) |>
    summarize(n_eea_10km=n(), my_total_area=sum(area))

##### range with custom function

my_species_range = list()
species_names <- unique(my_distribution_summary$Species)
source("necca_spatial_functions.R")
for (i in seq_along(species_names)) {
    dist <- my_distribution_whole_grid |>
        filter(Species==species_names[i])
    
    grids <- eea_10km_ETRS89 |>
        filter(CELLCODE %in% dist$CELLCODE) |>
        mutate(cell_origin="distribution") 

    my_species_range[[i]] <- expand_range_with_gap_distance(distribution = grids,
                                           full_grid = eea_10km_ETRS89,
                                           gap_distance_m = 40000, cellcode_col = "CELLCODE") |>
    mutate(Species=species_names[i])


}

my_species_range_all <- bind_rows(my_species_range)

## points over natura
points_n2k_joined <- st_join(points_u,
                             N2000_v32_ETRS89,
                             join = st_intersects,
                             left = TRUE,
                             largest = FALSE)

points_n2k <- st_intersects(points_n2k_joined, N2000_v32_ETRS89)
# do it at once with mutate
points_distinct_n2k <- points_n2k_joined |> 
    mutate(n2k = lengths(points_n2k) > 0)

points_distinct_n2k_only <- points_distinct_n2k |>
    filter(n2k==TRUE)

# my population over natura
my_population_n2k <- st_join(points_distinct_n2k_only,
                             eea_1km_ETRS89,
                             join = st_intersects,
                             left = TRUE,
                             largest = FALSE)


my_population_n2k_summary <- my_population_n2k |>
    st_drop_geometry() |> 
    distinct(Species,CELLCODE) |>
    group_by(Species) |>
    summarise(my_population_n2k=n(), .groups="keep")

################ all my calculations ################
summary_list <- list(my_distribution_summary, my_population_n2k_summary,my_population_summary)

my_summaries <- reduce(summary_list, ~ left_join(.x, .y, by = "Species"))

write_delim(my_summaries, "../anadoxos_deliverables/R_summaries.tsv",delim="\t")
#### summary
#### distribution
distribution_summary <- distribution |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    st_drop_geometry() |>
    group_by(SPECIES) |>
    summarise(n_eea_10km_cells=n(),
              total_dist_area=round(drop_units(sum(area)),digits=2))

distribution_unique <- distribution |>
    distinct(SPECIES,CELLCODE) |>
    group_by(SPECIES) |>
    summarise(unique_cells=n())

distribution_summary <- distribution_summary |>
    left_join(distribution_unique)

write_delim(distribution_summary, "../anadoxos_deliverables/validation_distribution.tsv",delim="\t")

print("check the distribution unique cells")
sum(distribution_summary$n_cells) == sum(distribution_summary$unique_cells)

print("check the distribution cells with the number of features")
sum(distribution_summary$n_cells)==nrow(distribution)

#### range summary
range_summary <- ranges |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    st_drop_geometry() |>
    group_by(SPECIES) |>
    summarise(total_range_area=round(drop_units(sum(area)),digits=2))

write_delim(range_summary, "../anadoxos_deliverables/validation_range.tsv",delim="\t")

#### check distinct cells
#### in range the cells of distribution must be present
range_dist_summary <- ranges |>
    distinct(SPECIES,CELLCODE,geometry) |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    st_drop_geometry() |>
    group_by(SPECIES) |>
    summarise(n_range_cells=n(),
              total_range_area=round(drop_units(sum(area)),digits=2))

range_dist_check <- ranges |>
    distinct(SPECIES,CELLCODE) |>
    group_by(SPECIES) |>
    summarise(unique_cells=n())

range_dist_summary <- range_dist_summary |>
    left_join(range_dist_check)

print("check the distribution unique cells")
sum(range_dist_summary$n_range_cells) == sum(range_dist_summary$unique_cells)

print("check the distribution cells with the number of features")
sum(range_dist_summary$n_range_cells)==nrow(ranges)

## check cells in range and dist
## join the files and calculate the number of files
## of each species AND CELLCODE
## the rows with 1 file should be the extra cells of 
## range.
## the rows with 2 files must be the rows of distribution
distribution_range_check <- rbind(distribution,ranges) |>
    dplyr::distinct(CELLCODE,SPECIES,file) |>
    group_by(CELLCODE,SPECIES) |>
    summarise(files=n_distinct(file),
              file_name=paste(file, collapse = ", "),
              .groups = "drop")


print("check that the distribution cells are in range")

distribution_range_joined <- distribution_range_check |>
    filter(files>1)
print("Are all cells that are in both files the same size as in distribution file?")
nrow(distribution_range_joined) == nrow(distribution)

print("Are all cells of combind in range file?")
sum(range_dist_summary$n_cells)==nrow(distribution_range_check)


#sf_use_s2(TRUE)
#### population, calculate 

population_summary <- population |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    mutate(area_tr=units::set_units(st_area(st_transform(geometry,4326)),"km^2")) |>
    st_drop_geometry() |>
    group_by(Species) |>
    summarise(
              total_1km_area=round(drop_units(sum(area)),digits=2),
              total_1km_area_tr=round(drop_units(sum(area_tr)),digits=2)) |>
    rename("SPECIES"="Species")


#### population check distinct population 
population_distinct <- population |>
    distinct(Species,CELLCODE_1,geometry)
print("check the population unique cells")
nrow(population) == nrow(population_distinct)

st_write(population_distinct,
         "../anadoxos_deliverables/population_distinct.gpkg",
         layer = "population",
         delete_layer = TRUE)
#### there are duplicate cells 

population_distinct_area <- population |>
    distinct(Species,CELLCODE_1,geometry) |>
    mutate(area=units::set_units(st_area(geometry),"km^2")) |>
    mutate(area_tr=units::set_units(st_area(st_transform(geometry,4326)),"km^2")) |>
    st_drop_geometry() |>
    group_by(Species) |>
    summarise(
              n_eea_1km_cells=n(),
              total_1km_u_area=round(drop_units(sum(area)),digits=2),
              total_1km_u_area_tr=round(drop_units(sum(area_tr)),digits=2)) |>
    left_join(population_summary, by=c("Species"="SPECIES"))

write_delim(population_distinct_area, "../anadoxos_deliverables/validation_population.tsv",delim="\t")

#### population over N2K

population_n2k_distinct <- population_n2k |>
    distinct(Species, CELLCODE_1)


print("check the population n2k unique cells")
nrow(population_n2k) == nrow(population_n2k_distinct)

population_n2k_summary <- population_n2k |>
    distinct(Species, CELLCODE_1) |>
    group_by(Species) |>
    summarise(population_n2k=n()) |>
    rename("SPECIES"="Species")


######################### Check my calculation with Deliverables ################
deliverable_summary_list <- list(range_summary, distribution_summary, population_n2k_summary,population_summary)

deliverable_summaries <- reduce(deliverable_summary_list, ~ left_join(.x, .y, by = "SPECIES"))

write_delim(deliverable_summaries, "../anadoxos_deliverables/deliverable_summaries_with_R.tsv",delim="\t")


######################################################### validate with map


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
