#!/usr/bin/Rscript

## Script name: species_enrichment.R
##
## Purpose of script: use public databases to retrieve information of 
## species regarding their taxonomy, their global distribution and their 
## IUCN status
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-11-06

library(sf)
library(tidyverse)
library(readxl)
library(units)
library(rnaturalearth)

############################# Load data ############################
## borders for maps
greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")

########################## species information ##############################
# the species of the previous epopteia 2015

invertebrates_92_43 <- readxl::read_excel("../data/invertebrates_92_43.xlsx", skip=2, col_names=F)

colnames(invertebrates_92_43) <- c("no","area","SPECIES_NAME","SPECIES_ID","ORDER","CLASS","PRIORITY","ANNEX_II","ANNEX_IV","ANNEX_V","KD","POPULATION_TREND","POPULATION_SIZE_UNIT","OCCURRENCE","SD")

species_names <- unique(invertebrates_92_43$SPECIES_NAME) 
##################### Natura2000 v32 version #############################

natura_v32_species <- read_delim("../data/Natura2000DB_V32_species.tsv",delim="\t")
natura_v32_site <- read_delim("../data/Natura2000DB_V32_site.tsv",delim="\t")
natura_v32_region <- read_delim("../data/Natura2000DB_V32_region.tsv",delim="\t")
natura_v32_other_species <- read_delim("../data/Natura2000DB_V32_other_species.tsv",delim="\t")

groups_df <- data.frame(SPECIES_GROUP=c("R","B","A","M","I","F","P"),
                    groups_names=c("Reptiles", "Birds", "Amphibians", "Mammals", "Invertebrates", "Fish", "Plants"))
site_types_df <- data.frame(SITE_TYPE=c("A","B","C"), SITE_TYPE_NAME=c("SPAs","SCIs_SACs", "SPAs_and_SCIs_SACs"))

### EU DEM Greece
#eu_dem_gr <- rast("../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_4326_gr.tif")

### EU DEM Greece slope
#eu_dem_slope <- rast("../spatial_data/EU_DEM_slope_gr/crop_eudem_slop_3035_europe.tif")

### World Clim , Bioclim Variables
###
#world_clim_directory <- "../spatial_data/world_clim_greece/"
#world_clim_files <- list.files(world_clim_directory, pattern = "\\.tif$", full.names = TRUE)

#world_clim_list <- lapply(world_clim_files, rast)


### Ecosystem Types of Europe 3.1 Terrestrial
#ecosystem_types_gr <- rast("/Users/talos/Documents/programming_projects/necca_epopteia/spatial_data/ecosystem_types_gr/crop_eea_r_3035_100_m_etm-terrestrial-c_2012_v3-1_r00.tif")

### Vegetation_map_Greece
## the encoding ENCODING=WINDOWS-1253 helped to see the greek characters
#vegetation_map <- sf::st_read("../spatial_data/Vegetation_map_Greece/D_xabxg_VPG_60-98_GEO.shp",
#                              options = "ENCODING=WINDOWS-1253")

#vegetation_map_wgs <- st_transform(vegetation_map,4326) |>
#    st_make_valid()

#veg_data <- data.frame(
#  A_VEG_TYPE = c("ΕΛΑ", "ΕΡΛ", "ΠΜΑ", "ΠΛΔ", "ΠΔΑ", "ΠΧΑ", "ΠΚΟ", "ΠΘΑ", "ΚΠΡ",
#                 "ΑΡΚ", "ΟΞΥ", "ΔΡΥ", "ΚΑΣ", "ΣΗΜ", "ΣΦΕ", "ΦΙΛ", "ΦΠΛ", "ΠΑΡ",
#                 "ΕΥΚ", "ΦΟΙ", "ΘΑΜ", "ΦΘΑ", "ΛΙΒ", "ΑΓΟ", "ΟΙΚ", "ΓΚΑ", "ΓΚΕ",
#                 "ΛΧΡ", "ΛΙΜ"),
#  A_VEG_NAME = c("Ελάτη", "Ερυθρελάτη", "Πεύκη μαύρη", "Πεύκη λευκόδερμη", "Πεύκη δασική",
#                 "Πεύκη χαλέπιος", "Πεύκη κουκουναριά", "Πεύκη θαλασσία", "Κυπαρίσσι",
#                 "Άρκευθος", "Οξυά", "Δρύς", "Καστανιά", "Σημύδα", "Σφένδαμος",
#                 "Φιλύρα", "Φυλλοβόλα πλατύφυλλα", "Παραποτάμια βλάστηση", "Ευκάλυπτος",
#                 "Φοίνικας", "Θάμνοι", "Φυλλοβόλοι θάμνοι", "Λιβάδια, αραιά ξυλ. βλάστηση",
#                 "Άγονα", "Οικισμοί", "Γεωργ. καλλιέργειες", "Γεωργ. καλλιέργειες εγκατ.",
#                 "Λοιπές χρήσεις", "Λίμνη"),
#  stringsAsFactors = FALSE
#)
#
#### EUNIS_Habitats_2018 for the Natura2000 areas
## the same as 01_Χαρτογράφηση χερσαίων Τ.Ο._EKXA
#
#EUNIS_Habitats <- sf::st_read("../spatial_data/EUNIS_Habitats_2018/Habitats_2018/Habitats.shp")
#
#EUNIS_Habitats_wgs <- st_transform(EUNIS_Habitats,4326) |>
#    st_make_valid()
#
########## HILDA Greece land use change ##########
#hilda_cat <- data.frame(hilda_id = c("11","22","33","44","55","66","77"),
#                        hilda_name=c("urban","cropland","pasture/rangeland",
#                                     "forest", "unmanaged grass/shrubland","sparse/no vegetation", "water"),
#                        hilda_hex=c("#000000","#AE6120","#98BA6A","#07A07D","#BE81A3","#999999", "#1370A1"))
#hilda_cat_v <- c("urban"="#000000",
#                 "cropland"="#AE6120",
#                 "pasture/rangeland"="#98BA6A",
#                 "forest"="#07A07D",
#                 "unmanaged grass/shrubland"="#BE81A3",
#                 "sparse/no vegetation"="#999999",
#                 "water"="#1370A1")
################################ Descriptives ###################################
##  
all_species <- unique(c(unique(natura_v32_species$SPECIES_NAME),unique(natura_v32_other_species$OTHER_SPECIES_NAME)))

natura_v32_species_sum <- natura_v32_species |> 
    distinct(SPECIES_NAME,SPECIES_GROUP,SITE_CODE) 

natura_v32_other_species_sum <- natura_v32_other_species |> 
    distinct(OTHER_SPECIES_NAME,OTHER_SPECIES_GROUP, SITE_CODE)

colnames(natura_v32_other_species_sum) <- colnames(natura_v32_species_sum)

natura_v32_all_species <- rbind(natura_v32_other_species_sum, natura_v32_species_sum) |>
    distinct()

### species
groups_species_summary <- natura_v32_all_species |>
    distinct(SPECIES_NAME, SPECIES_GROUP) |>
    group_by(SPECIES_GROUP) |> 
    summarise(n_species=n()) |>
    left_join(groups_df)

invertebrates_all <- natura_v32_all_species |>
    filter(SPECIES_GROUP=="I")

invertebrates_all_natura_summary <- invertebrates_all |>
    group_by(SITE_CODE) |>
    summarise(n_species=n()) |> 
    left_join(natura_v32_region)

invertebrates_all_natura_summary_region <- invertebrates_all_natura_summary |>
    group_by(REGION_NAME) |>
    summarise(n_species=sum(n_species), n_sites=n())

invertebrates_all_species <- invertebrates_all |>
    distinct(SPECIES_NAME)

### regions
## some sites are in multiple regions 
# > natura_v32_region |> distinct(REGION_CODE,SITE_CODE) |> group_by(SITE_CODE) |> summarise(n=n()) |> arrange(desc(n))

## overlap
## 16 not included, 23 included
invertebrates_not_in_natura <- species_names[which(!(species_names %in% natura_v32_species$SPECIES_NAME))]
invertebrates_not_in_other_natura <- species_names[which(!(species_names %in% natura_v32_other_species$OTHER_SPECIES_NAME))]


invertebrates_not_everywhere <- species_names[which(!(species_names %in% all_species))] 

species_names %in% unique(natura_v32_other_species$OTHER_SPECIES_NAME)

## Species from excel statistics

natura_v32_all_species_excel <- natura_v32_all_species |>
    filter(SPECIES_NAME %in% species_names) |>
    group_by(SPECIES_NAME) |>
    summarise(n_sites=n(),SITE_CODE=str_c(SITE_CODE, collapse = ","))

write_delim(natura_v32_all_species_excel, "../results/natura_v32_all_species_excel.tsv", delim="\t")

### sites

groups_sites_summary <- natura_v32_all_species |>
    distinct(SITE_CODE, SPECIES_GROUP,SPECIES_NAME) |>
    group_by(SITE_CODE,SPECIES_GROUP) |> 
    summarise(n_species=n(),
              species_names=str_c(SPECIES_NAME,collapse = ","),
              .groups="keep") |>
    left_join(natura_v32_site) |>
    left_join(site_types_df)

write_delim(groups_sites_summary,
            "../results/natura_groups_summary_all.tsv",
            delim="\t")

groups_sites_summary_invertebrates <- groups_sites_summary |>
    filter(SPECIES_GROUP=="I")

write_delim(groups_sites_summary_invertebrates,
            "../results/natura_groups_summary_invertebrates.tsv",
            delim="\t")

#### sites
groups_site_type_summary <- natura_v32_all_species |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME, SPECIES_GROUP,SITE_CODE) |>
    left_join(groups_df) |>
    group_by(SITE_TYPE_NAME,groups_names) |>
    summarise(n_sites=n(), .groups="keep") 

groups_site_type_sites_summary_all <- natura_v32_all_species |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME, SITE_CODE) |>
    group_by(SITE_TYPE_NAME) |>
    summarise(n_sites=n()) |>
    mutate(groups_names="All") |>
    rbind(groups_site_type_summary)

#### type species
groups_site_type_species_summary <- natura_v32_all_species |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME, SPECIES_GROUP,SPECIES_NAME) |>
    left_join(groups_df) |>
    group_by(SITE_TYPE_NAME,groups_names) |>
    summarise(n_species=n(), .groups="keep") 

groups_site_type_species_summary_all <- natura_v32_all_species |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME,SPECIES_NAME) |>
    group_by(SITE_TYPE_NAME) |>
    summarise(n_species=n()) |>
    mutate(groups_names="All") |>
    rbind(groups_site_type_species_summary)


### species from annex II
groups_sites_invertebrates <- groups_sites_summary |>
    filter(SPECIES_GROUP=="I")

species_annex_II <- invertebrates_92_43 |>
    filter(!is.na(ANNEX_II))

groups_site_type_species_annex_II <- natura_v32_all_species |>
    filter(SPECIES_NAME %in% species_annex_II$SPECIES_NAME) |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME, SPECIES_GROUP,SPECIES_NAME) |>
    mutate(groups_names="Invertebrates_annex_II") |>
    group_by(SITE_TYPE_NAME,groups_names) |>
    summarise(n_species=n(), .groups="keep") 

groups_site_type_annex_II_summary <- natura_v32_all_species |>
    filter(SPECIES_NAME %in% species_annex_II$SPECIES_NAME) |>
    left_join(natura_v32_site) |>
    left_join(site_types_df) |>
    distinct(SITE_TYPE_NAME, SPECIES_GROUP,SITE_CODE) |>
    mutate(groups_names="Invertebrates_annex_II") |>
    group_by(SITE_TYPE_NAME,groups_names) |>
    summarise(n_sites=n(), .groups="keep") 

natura_sites_groups_species_annex_II <- groups_site_type_species_annex_II |>
    left_join(groups_site_type_annex_II_summary,
    by=c("SITE_TYPE_NAME","groups_names"))

### sites and species together

natura_sites_groups_species_summary <- groups_site_type_sites_summary_all |>
    left_join(groups_site_type_species_summary_all,
              by=c("SITE_TYPE_NAME","groups_names")) |>
    rbind(natura_sites_groups_species_annex_II)

write_delim(natura_sites_groups_species_summary,
            "../results/natura_sites_groups_species_summary.tsv",
            delim="\t")


### rasters stack
#rasters_list <- list(eu_dem_gr,eu_dem_slope,ecosystem_types_gr) # hilda_rast_list, add later ,world_clim_list

#results_ext <- extract_from_named_rasters(rasters_list,points_sf)
#results_ext_2 <- extract_from_named_rasters(world_clim_list,results_ext)

### shapefiles

####################### extract ####################### 
#result <- extract_polygon_info_multi(
#  results_ext_2,
#  polygons_list <- list(eea_1km_wgs,
#                      eea_10km_wgs,
#                      N2000_v32_wgs_sci,
#                      N2000_v32_wgs_spa,
#                      N2000_v32_wgs_scispa,
#                      vegetation_map_wgs,
#                      EUNIS_Habitats_wgs),
#  suffixes = c("eea_1km","eea_10km","N2000_v32_sci","N2000_v32_spa","N2000_v32_scispa","vegetation_map","EUNIS_Habitats")
#)
#
#species_occurrences_spatial <- result
#write_delim(species_occurrences_spatial,"../results/species_occurrences_spatial.tsv",delim="\t")
#
############################## Hilda analysis ##############################

##### Hilda land cover change
#hilda_path <- "../spatial_data/hilda_greece/"
#hilda_id_names <- read_delim(paste0(hilda_path, "hilda_transitions_names.tsv", sep=""), delim="\t")
#hilda_all <- list.files(hilda_path)
#hilda_files <- hilda_all[grepl("*.tif", hilda_all)] 
#hilda_rast_list <- lapply(hilda_files, function(x) rast(paste0(hilda_path,x,sep="")))
#
results_ext_3 <- extract_from_named_rasters(hilda_rast_list,points_final)


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
