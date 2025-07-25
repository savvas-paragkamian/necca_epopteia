#!/usr/bin/Rscript

## Script name: data_validation_qc.R
##
## Purpose of script:
## In this script we perform validation and quality control
## of the deliverables of the project. At the species level 
## the population (number of 1X1 eea cells), distribution and n2k dist.
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
species_samples_art17 <- st_read("../results/species_samples_art17_all.tsv")

wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")

greece_regions_ETRS89 <- st_transform(greece_regions,3035)

eea_1km <- sf::st_read("../spatial_data/eea_reference_grid/gr_1km.shp")
eea_1km_ETRS89 <- st_transform(eea_1km, 3035)

eea_10km <- sf::st_read("../spatial_data/eea_reference_grid/gr_10km.shp")
eea_10km_ETRS89 <- st_transform(eea_10km, 3035)
eea_10km_wgs <- st_transform(eea_10km,4326)

## do the st_union because greece_regions_ETRS89 is a multipolygon!!
eea_10km_ETRS89_gr <- st_intersection(eea_10km_ETRS89, st_union(greece_regions_ETRS89)) |>
    distinct(CELLCODE,geometry) |> 
    mutate(area=units::set_units(st_area(geometry),"km^2"))

eea_10km_ETRS89_gr_area <- st_drop_geometry(eea_10km_ETRS89_gr)

N2000_v32 <- sf::st_read("../spatial_data/N2000_spatial_GR_2021_12_09_v32/N2000_spatial_GR_2021_12_09_v32.shp")

N2000_v32_wgs <- st_transform(N2000_v32,4326)
N2000_v32_ETRS89 <- st_transform(N2000_v32, 3035)

N2000_v32_ETRS89_sci <- N2000_v32_ETRS89 |>
    filter(SITETYPE!="SPA")

## elevation
### EU DEM Greece
eu_dem_gr <- rast("../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_4326_gr.tif")

### EU DEM Greece slope
eu_dem_slope <- rast("../spatial_data/EU_DEM_slope_gr/crop_eudem_slop_3035_europe.tif")


########################### Validation - QC ###########################
######### Anadohos ExtendedTemplate_Reporting dataSpecies_total ##################
#### latest check 2025-07-25
#### load files

points <- sf::st_read("../anadoxos_deliverables/Maps/Valid species/VerifiedOccurrenceDB_PlusOrphans_LAEA.shp")

ranges <- sf::st_read("../anadoxos_deliverables/Maps/Valid species/Species_range_invertebrates.shp") |>
        mutate(file="range")

distribution <- sf::st_read("../anadoxos_deliverables/Maps/Valid species/Species_distribution_invertebrates.shp") |>
    mutate(file="distribution")

population <- sf::st_read("../anadoxos_deliverables/Maps/Valid species/National_gr1km_Pops.shp")

population_n2k <- sf::st_read("../anadoxos_deliverables/Maps/Valid species/N2K_gr1km_Pops.shp")

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

my_population_unio <- my_population |>
    filter(grepl("Unio*.",Species)) |>
    filter(Species!="Unio elongatulus") |>
    mutate(Species_complex = "Unio crassus complex") |>
    distinct(Species_complex,CELLCODE) |>
    group_by(Species_complex) |>
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

##### troubleshooting the range function! not the same as in range tool
#
#    distribution <- my_distribution_whole_grid |>
#        filter(Species=="Apatura metis")
#    grids <- distribution # species distribution
#
#    # 1. Calculate centroid distances
#    centroids <- st_centroid(grids)
#    dist_matrix <- st_distance(centroids, which="Euclidean")
#    dist_num <- drop_units(dist_matrix)
#    
#    # 2. Name rows/columns of matrix by CELLCODE
#    cellcode <- grids[["CELLCODE"]]
#    rownames(dist_num) <- colnames(dist_num) <- cellcode
#    
#    # 3. Create adjacency matrix: keep only upper triangle and filter by gap
#    upper_only <- dist_num
#    upper_only[!upper.tri(upper_only)] <- 0
#    adj_matrix <- upper_only
#    adj_matrix[adj_matrix <= 0 | adj_matrix > 400000] <- 0
#    diag(adj_matrix) <- 0
#    
#    # 4. Convert adjacency matrix to edgelist
#    edge_list <- which(adj_matrix != 0, arr.ind = TRUE)
#    edge_df <- data.frame(
#        from = rownames(adj_matrix)[edge_list[, 1]],
#        to   = colnames(adj_matrix)[edge_list[, 2]],
#        weight = adj_matrix[edge_list]
#    )
#    
#    # 5. Join centroid geometries for from/to
#    centroids_df <- centroids %>% st_drop_geometry() %>% mutate(geometry = st_geometry(centroids))
#    from_geom <- centroids[match(edge_df$from, cellcode), ]
#    to_geom   <- centroids[match(edge_df$to, cellcode), ]
#    
#    # 6. Create lines between centroids
#    line_list <- mapply(
#      function(p1, p2) st_linestring(matrix(c(st_coordinates(p1), st_coordinates(p2)), ncol = 2, byrow = TRUE)),
#      st_geometry(from_geom),
#      st_geometry(to_geom),
#      SIMPLIFY = FALSE
#    )
#    
#    lines_sf <- st_sf(
#        from = edge_df$from,
#        to   = edge_df$to,
#        geometry = st_sfc(line_list, crs = st_crs(grids))
#    )
#    
#    # 7. Find cells in full grid intersecting these lines
#    full_grid <- eea_10km_ETRS89
#    intersecting_cells <- st_intersects(full_grid, lines_sf, sparse = FALSE)
#    filled_cells <- full_grid[apply(intersecting_cells, 1, any), ] %>%
#        mutate(cell_origin = "range")
#    
#    # 8. Return union of filled cells and original grids
#    expanded_grid <- rbind(filled_cells, grids) %>% distinct()
#    
#
#####
my_species_range_all <- bind_rows(my_species_range)

my_species_range_summary <- my_species_range_all |>
    st_drop_geometry() |> 
    distinct(Species,CELLCODE) |>
    left_join(eea_10km_ETRS89_gr_area) |>
    group_by(Species) |>
    summarise(my_eea10km_range=n(), my_range_area=sum(area,na.rm = TRUE))

########## points over natura
points_u_v <- points |> 
    filter(Distributi=="Verified") |>
#    st_drop_geometry() |>
    count(Species, decimalLat,decimalLon, geometry)


points_n2k_joined <- st_join(points_u_v,
                             N2000_v32_ETRS89_sci,
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
    filter(!is.na(SITETYPE)) |>
    distinct(Species,CELLCODE) |>
    group_by(Species) |>
    summarise(my_population_n2k=n(), .groups="keep")

my_population_n2k_unio <- my_population_n2k |>
    filter(grepl("Unio*.",Species)) |>
    filter(Species!="Unio elongatulus") |>
    mutate(Species_complex = "Unio crassus complex") |>
    distinct(Species_complex,CELLCODE) |>
    group_by(Species_complex) |>
    summarize(my_n_eea_1km=n())
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

########################### Validation - QC 2.6 Habitats ##################

habitats_2_6 <- sf::st_read("../anadoxos_deliverables/2.6. Αποτύπωση των ενδιαιτημάτων/Species_habitats_invertebrates/Species_habitats_invertebrates.shp")
habitats_prefs <- read_excel("../anadoxos_deliverables/2.6. Αποτύπωση των ενδιαιτημάτων/Παράρτημα Ι.xls")

points <- sf::st_read("../anadoxos_deliverables/Distribution&RangeMaps_20250611/VerifiedOccurrenceDB_PlusOrphans_LAEA.shp")

points_sf <- st_transform(points, crs(eu_dem_gr))

points_dem <- terra::extract(eu_dem_gr, vect(points_sf))

points_with_vals <- cbind(points_sf, points_dem[,-1]) 

points_dem_sum <- points_with_vals |>
    group_by(Species) |>
    summarise(dem_mean=mean(points_dem....1.,na.rm=T)) |>
    left_join(habitats_prefs, by=c("Species"="SPECIES"))




########################### Validation - QC Protogeni epopteia 2025 ##################

######## Deigmata 
deigmata_data <- read_xlsx("../anadoxos_deliverables/FINAL Invertebrates ΠΒΔ_V6.xlsx",
                           sheet="Δείγματα Ασπόνδυλων",
                           col_names=T
                           ) |> slice(-1)

deigmata_data_qc <- deigmata_data |> 
    filter(!is.na(Sam_ID)) |>
    mutate(latitude=as.numeric(`Γεωγραφικό Πλάτος (WGS84) Αρχη`),
           longitude=as.numeric(`Γεωγραφικό Μήκος (WGS84) Αρχή`)) 

########### qc missing coordinates

deigmata_data_qc_missing <- deigmata_data_qc |>
    filter(is.na(latitude) | is.na(longitude))


deigmata_data_sf <- deigmata_data_qc |>
    filter(!is.na(longitude)) |> 
    st_as_sf(
             coords=c("longitude","latitude"),
             remove=F,
             crs="WGS84")

st_write(deigmata_data_sf,
         "../anadoxos_deliverables/anadoxos_samples/deigmata_data.gpkg",
         layer = "deigmata_data_sf",
         delete_layer = TRUE)


######## Eidi

eidi_data <- read_xlsx("../anadoxos_deliverables/FINAL Invertebrates ΠΒΔ_V6.xlsx",
                       sheet="Είδη",
                       col_names=T
                           ) |> slice(-1) |> filter(!is.na(Obs_ID))

deigmata_eidi <- eidi_data |>
    mutate(species=if_else(`Όνομα είδους`=="Άλλο",
                           `Άλλο είδος`,
                           `Όνομα είδους`)) |>
    mutate(art17_92_43_EEC=if_else(`Όνομα είδους`!="Άλλο",
                                   TRUE,
                                   FALSE)) |>
    left_join(deigmata_data_qc, by=c("Sam_ID"="Sam_ID"))

#### export to gpkg
deigmata_eidi_sf <- deigmata_eidi |>
    filter(!is.na(longitude)) |> 
    st_as_sf(
             coords=c("longitude","latitude"),
             remove=F,
             crs="WGS84")

st_write(deigmata_eidi_sf,
         "../anadoxos_deliverables/anadoxos_samples/deigmata_eidi.gpkg",
         layer = "deigmata_eidi_sf",
         delete_layer = TRUE)

##### species from the list
deigmata_eidi_art17 <- deigmata_eidi |>
    distinct(Obs_ID,Sam_ID, species, art17_92_43_EEC) |>
    filter(art17_92_43_EEC==TRUE)

sort(unique(deigmata_eidi_art17$species))

deigmata_eidi_art17_summary <- deigmata_eidi_art17 |>
    group_by(species,Sam_ID) |>
    summarise(n_obs=n(),.groups="keep")
    
print("number of samples that have species from art 17")
length(unique(deigmata_eidi_art17_summary$Sam_ID))

print("number of observations that have species from art 17")
sum(deigmata_eidi_art17_summary$n_obs)

deigmata_eidi_obs <- deigmata_eidi |>
    group_by(species,Sam_ID) |>
    summarise(n_obs=n(),.groups="keep")

deigmata_eidi_obs_s <- deigmata_eidi_obs |>
    group_by(species) |>
    summarise(n_obs=sum(n_obs),n_samples=n())


# which samples don't have occurrences
samples_no_species <- deigmata_data_qc |> 
    filter(!(Sam_ID %in% unique(deigmata_eidi$Sam_ID)))

samples_with_species <- deigmata_data_qc |> 
    mutate(sample_with_species = Sam_ID %in% unique(deigmata_eidi$Sam_ID)) |>
    mutate(samples_with_art17 = Sam_ID %in% unique(deigmata_eidi_art17_summary$Sam_ID)) |>
    mutate(species_category = if_else(sample_with_species==FALSE,"without species",
                                      if_else(samples_with_art17==TRUE,
                                              "species from list",
                                              "other species")))

print("summary of samples with species")
table(samples_with_species$sample_with_species)
print("summary of samples with species from list")
table(samples_with_species$samples_with_art17)
print("summary of samples category based on species")
table(samples_with_species$species_category)


################### eea 10km grid overlap ################

deigmata_data_sf_ETRS89 <- st_transform(deigmata_data_sf,3035)
points_eea_joined <- st_join(deigmata_data_sf_ETRS89,
                             eea_10km_ETRS89,
                             join = st_intersects,
                             left = TRUE,
                             largest = FALSE)

points_eea_joined_distinct <- points_eea_joined |>
    distinct(Sam_ID, CELLCODE) |>
    left_join(eea_10km_wgs) |>
    sf::st_as_sf()

print("number of eea 10 sp km grids")
length(unique(points_eea_joined_distinct$CELLCODE))

#### map

epopteia_2025_samples_with_species_map <- ggplot() +
    geom_sf(greece_regions, mapping=aes()) +
    geom_sf(points_eea_joined_distinct, mapping=aes(alpha=0.8),color="lightseagreen" ) +
    geom_point(samples_with_species,
            mapping=aes(x=longitude,
                        y=latitude,
                        color=species_category),
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

ggsave("../anadoxos_deliverables/anadoxos_samples/epopteia_2025_samples_with_species_map.png", 
       plot=epopteia_2025_samples_with_species_map, 
       height = 20, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="png")


################### greece overlap ################
dist_matrix <- st_distance(deigmata_data_sf, greece_regions)
deigmata_data_sf$min_dist_m <- apply(dist_matrix, 1, min)

# points inside
points_inside_or_touching <- deigmata_data_sf[deigmata_data_sf$min_dist_m ==0, ]

points_outside <- deigmata_data_sf[deigmata_data_sf$min_dist_m >0,]

print("Samples with coordinates in the sea or out of hellenic border")
deigmata_data_sf |> filter(min_dist_m>0) |> dplyr::select(Sam_ID, min_dist_m) 

epopteia_2025_species_gr_map <- ggplot() +
    geom_sf(greece_regions, mapping=aes()) +
    geom_point(deigmata_data_sf,
            mapping=aes(x=longitude,
                        y=latitude,
                        color=`Χρηματοδοτικό Μέσο`),
            size=1.8,
            alpha=0.8,
            show.legend=T) +
    geom_point(points_outside,
            mapping=aes(x=longitude,
                        y=latitude,
                        color="outside"),
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

ggsave("../anadoxos_deliverables/anadoxos_samples/epopteia_2025_species_gr_map.png", 
       plot=epopteia_2025_species_gr_map, 
       height = 20, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="png")

################### n2k overlap ################

points_n2k_joined <- st_join(deigmata_data_sf_ETRS89,
                             N2000_v32_ETRS89,
                             join = st_intersects,
                             left = TRUE,
                             largest = FALSE)

points_n2k <- st_intersects(points_n2k_joined, N2000_v32_ETRS89)

# do it at once with mutate
points_distinct_n2k <- points_n2k_joined |> 
    mutate(n2k = lengths(points_n2k) > 0) |>
    distinct(Sam_ID,n2k,longitude,latitude) |>
    left_join(samples_with_species)

print("species category per samples and appearrance in NATURA2000")
table(points_distinct_n2k$n2k,points_distinct_n2k$species_category)


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
### natura2000 with all points
g_art17_n2000 <- g_base_n2000 +
    geom_point(points_distinct_n2k,
            mapping=aes(x=longitude,
                        y=latitude,
                        color=n2k),
            size=1.2,
            alpha=0.8,
            show.legend=T) +
#    scale_color_manual(values=datasets_colors,
#                        name = "Datasets")+
    guides(
           fill=guide_legend(position = "inside",override.aes = list(linetype = 0,color=NA)),
           color=guide_legend(position = "inside",override.aes = list(linetype = 0,fill=NA)))+
    theme(legend.position.inside = c(0.87, 0.75)
    )


ggsave("../anadoxos_deliverables/anadoxos_samples/samples_invertebrates_natura.png", 
           plot=g_art17_n2000, 
           height = 20, 
           width = 25,
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
