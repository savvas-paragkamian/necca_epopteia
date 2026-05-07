library(sf)
library(terra)
library(units)
library(readxl)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggnewscale)



#----------------------------------------------------------------------------#
################################# MAPS #######################################
##############################################################################
#----------------------------------------------------------------------------#


# plot

points_gr_plot <- ggplot() +
    geom_sf(data = eea_10km_ETRS89, aes(color = "eea"), fill="transparent", shape = 16, show.legend = TRUE) +
    geom_sf(data = greece_regions_ETRS89, aes(color = "Polygons"), fill = NA, show.legend = TRUE) +
    geom_sf(data = points_inside_or_touching, aes(color = "Points inside"), shape = 1, show.legend = TRUE) +
    geom_sf(data = points_500m, aes(color = "Points away < 500m"), shape = 1, show.legend = TRUE) +
    geom_sf(data = points_away_500m, aes(color = "Points away > 500m"), shape = 1, show.legend = TRUE) +
    scale_color_manual(
      name = "Feature Type",
      values = c(
        "Gr borders" = "blue",
        "Points inside" = "black",
        "Points away > 500m" = "blue",
        "Points away < 500m" = "red"
      )
    ) +
    theme_bw()+
    theme(legend.position="inside",
          legend.position.inside = c(0.87, 0.75))

ggsave(plot=points_gr_plot,
       "../figures/map_occurrences_borders_filtering.png",
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")





species_with_data <- unique(species_samples_presence_final_sf$species)
datasets_colors <- c(
                     "GBIF"="seagreen",
                     "NECCA_redlist"="#B31319",
                     "E1X_MDPP_2014_2024"="#FDF79C",
                     "E1X_DB"="#2BA09F",
                     "E2X_DB"="maroon1",
                     "db_refs_necca_2025"="green1",
                     "E1X_DB_references"="#141D43",
                     "DistrMap_2013_2018"="darkorchid1",
                     "Action Plan 2019"="#F85C29"
                     )

# base plot with Natura2000 areas of Greece
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
    geom_sf(species_samples_presence_final_sf,
            mapping=aes(
                        shape=basisOfRecord,
                        color=datasetName),
            size=2,
            alpha=0.8,
            show.legend=T) +
    scale_color_manual(values=datasets_colors,
                        name = "Datasets")+
    guides(
           fill=guide_legend(position = "right",override.aes = list(linetype = 0,color=NA)),
           shape=guide_legend(position = "right"),
           color=guide_legend(position = "right",override.aes = list(linetype = 0,fill=NA)))+
    theme(legend.position = "right")

ggsave("../figures/map_art17_invertebrates_natura.png", 
           plot=g_art17_n2000, 
           height = 30, 
           width = 40,
           dpi = 300, 
           units="cm",
           device="png")

#----------------------------------------------------------------------------#
############################ Species Occurrences #############################
#----------------------------------------------------------------------------#

### figures of each invertebrate of art17 for Greece

for (i in seq_along(species_with_data)){
    species_occurrences <- species_samples_presence_final_sf |>
        filter(species==species_with_data[i])

    dataset_colors_f <- datasets_colors[unique(species_occurrences$datasetName)]
    print(species_with_data[i])

    
    species_gr_map <- g_base_n2000 +
        geom_sf(species_occurrences,
                mapping=aes(
                        shape=basisOfRecord,
                        color=datasetName),
                size=2,
                alpha=0.9,
                show.legend=T) +
        scale_color_manual(values=dataset_colors_f,
                        name = "Datasets") +
        guides(
               fill=guide_legend(position = "right",override.aes = list(linetype = 0,color=NA)),
               shape=guide_legend(position = "right"),
               color=guide_legend(position = "right",override.aes = list(linetype = 0,fill=NA)))+
        ggtitle(paste(species_with_data[i]))+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position = "right",
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_occurrences.png", sep=""), 
           plot=species_gr_map, 
           height = 28, 
           width = 35,
           dpi = 300, 
           units="cm",
           device="png")

}

#----------------------------------------------------------------------------#
################################ Species Range ##############################
#----------------------------------------------------------------------------#
# calculate range for all species

species_range = list()
polygons_no_points_all <- list()
colors_cell_origin <- c("range"="mediumvioletred",
                        "distribution"="limegreen")

for (i in seq_along(species_with_data)){
    #i=1
    species_occurrences <- species_samples_presence_final_sf |>
        filter(species==species_with_data[i])
    species_dist <- species_samples_presence_final_sf |>
        filter(species==species_with_data[i]) |>
        filter(includeDistribution==TRUE) |>
        st_drop_geometry() |>
        distinct(species, CELLCODE_eea_10km) |>
        rename("CELLCODE"="CELLCODE_eea_10km")
    
    grids <- eea_10km_ETRS89 |>
        filter(CELLCODE %in% unique(species_dist$CELLCODE))
    
    # run the function for range calculations
    expanded_range <- expand_range_with_gap_distance(distribution = grids,
                                                     full_grid = eea_10km_ETRS89,
                                                     gap_distance_m = 40000)

    species_range[[species_with_data[i]]] <- expanded_range

    cell_origin_colors_f <- colors_cell_origin[unique(expanded_range$cell_origin)]

    # Maps
    species_gr_map <- g_base_n2000 +
        new_scale_fill()+
        geom_sf(expanded_range, mapping=aes(fill=cell_origin),
                alpha=0.8,
                colour="transparent",
                na.rm = F) +
        scale_fill_manual(
                          values=cell_origin_colors_f ,
                          guide = guide_legend(
                                                override.aes = list(
                                                                    linetype="solid",
                                                                    shape = NA)
                                                ),
                           name="Cell origin"
                           )+
        geom_sf(species_occurrences,
                mapping=aes(
                            color=datasetName),
                size=1,
                alpha=0.9,
                show.legend=T) +
        scale_color_manual(values=dataset_colors_f,
                            name = "Datasets")+
        guides(
               fill=guide_legend(position = "right",override.aes = list(linetype = 0,color=NA)))+
        ggtitle(paste(species_with_data[i]))+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position = "right",
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_range.png", sep=""), 
           plot=species_gr_map, 
           height = 28, 
           width = 35,
           dpi = 300, 
           units="cm",
           device="png")
    
}

all_range <- bind_rows(
  lapply(names(species_range), function(name) {
    species_range[[name]] %>%
      mutate(species = name)
  })
)

st_write(all_range, "../results/species_range/species_range.shp",append=TRUE)

##### calculate range

all_range_species <- all_range |>
    st_drop_geometry() |>
    distinct(CELLCODE,species) |>
    group_by(species) |>
    summarise(n_cells=n())

#----------------------------------------------------------------------------#
############################## Species Distribution ##############################
#----------------------------------------------------------------------------#

for (i in seq_along(species_with_data)){

    species_occurrences <- species_samples_presence_final_sf |>
        filter(species==species_with_data[i])

    species_dist <- species_samples_presence_final_sf |>
        filter(species==species_with_data[i]) |>
        filter(includeDistribution==TRUE) |>
        st_drop_geometry() |>
        distinct(species, CELLCODE_eea_10km) |>
        rename("CELLCODE"="CELLCODE_eea_10km")
    
    # summary of points with grids
    locations_10_grid_samples <- species_occurrences |>
        distinct(CELLCODE_eea_10km, decimalLatitude, decimalLongitude) |>
        group_by(CELLCODE_eea_10km) |>
        summarise(n_samples=n(),.groups="keep") |>
        rename("CELLCODE"="CELLCODE_eea_10km") |>
        left_join(eea_10km_ETRS89) |>
        st_as_sf()

    # colors of datasets present
    dataset_colors_f <- datasets_colors[unique(species_occurrences$datasetName)]

    print(species_with_data[i])

    # Maps
    species_gr_map <- g_base_n2000 +
        new_scale_fill()+
        geom_sf(locations_10_grid_samples, mapping=aes(fill=n_samples),
                alpha=0.5,
                colour="transparent",
                na.rm = F) +
        scale_fill_viridis_c(option = "viridis",
                             direction=1) +
        geom_sf(species_occurrences,
                mapping=aes(
                            color=datasetName),
                size=0.5,
                alpha=0.5) +
        scale_color_manual(values=dataset_colors_f,
                            name = "Datasets")+
        guides(
               color=guide_legend(position = "right",override.aes = list(linetype = 0,fill=NA)),
               fill = guide_colourbar(position = "right",
                                      alpha = 1,
                                      ticks = F,
                                      label = T,
                                      title="n_samples",
                                      title.vjust = 0.8,
                                      order = 1))+
        ggtitle(paste(species_with_data[i]))+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position = "right",
              legend.box.background = element_blank())
    
    ggsave(paste0("../figures/species_maps/map_", species_with_data[i], "_distribution.png", sep=""), 
           plot=species_gr_map, 
           height = 28, 
           width = 35,
           dpi = 300, 
           units="cm",
           device="png")

}
