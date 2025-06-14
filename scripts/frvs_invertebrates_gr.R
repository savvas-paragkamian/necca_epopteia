#!/usr/bin/env Rscript

## Script name: frvs_invertebrates_gr.R
##
## Purpose of script: Calculate the Favourite Reference Values
## of the invertebrates in Greece
##
## Author: Savvas Paragkamian
##
## Date Created: 2025-03-13

library(sf)
library(terra)
library(tidyverse)
library(ggnewscale)
library(readxl)
library(partition)
library(dismo)
library(units)

############################# Load Spatial Data ########################
wgs84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp")
#hellenic_borders_shp <- sf::st_read("../spatial_data/hellenic_borders/hellenic_borders.shp")

gr_1km <- sf::st_read("../spatial_data/eea_reference_grid/gr_1km.shp") |>
    st_transform(., crs="WGS84")

gr_10km <- sf::st_read("../spatial_data/eea_reference_grid/gr_10km.shp") |>
    st_transform(., crs="WGS84")

N2000_v32 <- sf::st_read("../spatial_data/N2000_spatial_GR_2021_12_09_v32/N2000_spatial_GR_2021_12_09_v32.shp")

N2000_v32_wgs <- st_transform(N2000_v32,4326)

########################### Load Species Data ###########################
### Species occurrences enriched ######
species_samples_art17 <- read_delim("../results/species_samples_art17.tsv", delim="\t")

species_taxonomy <- read_delim("../results/species_gbif_taxonomy_curated.tsv",delim="\t")
iucn_art17_invert_all <- read_delim("../results/iucn_art17_invert_all.tsv", delim="\t")
iucn_art17_invert_no_tax <- iucn_art17_invert_all |>
    dplyr::select(-ends_with("Name"),submittedName)

sspecies_dist_national_rep <- sf::st_read("../spatial_data/National report_2013_2018_shp/GR_Art17_species_distribution.shp")
species_dist_national_rep_sens <- sf::st_read("../spatial_data/National report_2013_2018_shp/GR_Art17_species_distribution_sensitive.shp")


########################## Flowchart for FRVs ##########################
# keep only one species name from the synonyms
species_info <- species_samples_art17 |>
    distinct(submittedName) |> 
    left_join(species_taxonomy, by=c("submittedName"="verbatim_name")) |>
    left_join(iucn_art17_invert_no_tax, by=c("submittedName"="submittedName"))
# species spatial 
species_art17_spatial <- species_samples_art17 |>
    left_join(species_info) |>
    filter(!is.na(decimalLongitude)) |> 
    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),
             remove=F,
             crs="WGS84")

########################## Parnassius apollo ###########################

parnassius_dist <- sf::st_read("../data/Parnassius apollo AP 2019/AP_Papollo_Distribution_LAEA.shp") |>
    st_transform(crs="WGS84")

p_apollo_points <- species_art17_spatial |>
    filter(species=="Parnassius apollo") |>
    mutate(A_VEG_TYPE_vegetation_map = as.factor(A_VEG_TYPE_vegetation_map))

apollo_mean <- p_apollo_points |>
    summarise(
    across(
           where(is.numeric),
           \(x) mean(x, na.rm = TRUE)
           )
    )

write_delim(apollo_mean,"../results/apollo_mean.tsv", delim="\t")

### hotspot
locations_10_grid_samples <- st_join(gr_10km, p_apollo_points, left=F) |>
    distinct(geometry,CELLCODE, decimalLatitude, decimalLongitude) |>
    group_by(geometry,CELLCODE) |>
    summarise(n_samples=n(),.groups="keep")


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

hotspot_apollo_map <- ggplot()+
    geom_sf(greece_regions, mapping=aes()) +
    geom_sf(N2000_v32_wgs, mapping=aes(fill=SITETYPE),
            alpha=0.8,
            #colour="transparent",
            na.rm = F,
            show.legend=T) +
    scale_fill_manual(
                      values= natura_colors,
                       guide = guide_legend(position = "inside",
                                            override.aes = list(alpha = 0.5,
                                                                linetype="solid",
                                                                color = NA)),
                       name="Natura2000")+
    new_scale_fill()+
    geom_sf(locations_10_grid_samples, mapping=aes(fill=n_samples),
            alpha=0.8,
            colour="transparent",
            na.rm = F,
            show.legend=T) +
    geom_sf(p_apollo_points, mapping=aes(color=datasetName),size=1,alpha=0.6) +
    scale_fill_gradient(low="gray50",
                        high="gray5",
                        guide = "colourbar")+
    scale_color_manual(values=datasets_colors,
                        name = "Datasets")+
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
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position.inside = c(0.87,0.65),
          legend.box.background = element_blank())

ggsave("../figures/hotspots_parnassius_apollo_map.png", 
       plot=hotspot_apollo_map, 
       height = 20, 
       width = 25,
       dpi = 300, 
       units="cm",
       device="png")

############## Modeling FRV ##################
apollo <- p_apollo_points |>
    st_drop_geometry() |>
    dplyr::select(decimalLatitude,decimalLongitude)


apollo_points <- vect(apollo, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

## pseudo-absences

## colinearity of raster data
library(corrr)
species_samples_art17_parnasious_num <- p_apollo_points |>
    st_drop_geometry() |>
    dplyr::select(where(is.numeric)) |>
    correlate()

cor_ff <- species_samples_art17_parnasious_num |>
    pivot_longer(-term, names_to="to_term", values_to="pearson") |>
    filter(abs(pearson) > 0.5) |>
    filter(term!=to_term) |>
    filter(if_all(where(is.character), ~ str_detect(., "X_wc2.1|X_eudem")))



############ environmental data
## raster paths
gr_1km_terra <- file.path("../spatial_data/eea_reference_grid/gr_1km_terra.tif")
slope <- file.path("../spatial_data/EU_DEM_slope_gr/crop_eudem_slop_3035_europe.tif")
dem <- file.path("../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_4326_gr.tif")
bio12 <- file.path("../spatial_data/world_clim_greece/crop_wc2.1_30s_bio_12.tif")
bio15 <- file.path("../spatial_data/world_clim_greece/crop_wc2.1_30s_bio_15.tif")
vegetation_map <- file.path("../spatial_data/Vegetation_map_Greece/D_xabxg_VPG_60-98_GEO.tif")

## load rasters
gr_1km_rast <- rast(gr_1km_terra)
slope_rast <- rast(slope)
dem_rast <- rast(dem)
bio12_rast <- rast(bio12)
bio15_rast <- rast(bio15)
veg_rast <- rast(vegetation_map)

# crop and resample 
dem_rast_c <- crop(dem_rast,ext(gr_1km_rast))
slope_rast_c <- crop(slope_rast,ext(gr_1km_rast))
bio12_rast_c <- crop(bio12_rast,ext(gr_1km_rast))
bio15_rast_c <- crop(bio15_rast,ext(gr_1km_rast))
veg_rast_c <- crop(veg_rast,ext(gr_1km_rast))


dem_rast_r <- resample(dem_rast_c,gr_1km_rast)
slope_rast_r <- resample(slope_rast_c,gr_1km_rast)
bio12_rast_r <- resample(bio12_rast_c,gr_1km_rast)
bio15_rast_r <- resample(bio15_rast_c,gr_1km_rast)
veg_rast_r <- resample(veg_rast_c,gr_1km_rast)


# Combine them into a single stack
env_stack <- c(dem_rast_r,slope_rast_r,bio12_rast_r,bio15_rast_r,veg_rast_r)
names(env_stack) <- c("eudem_dem_4258_europe","eudem_slop_3035_europe","wc2.1_30s_bio_12", "wc2.1_30s_bio_15", "A_VEG_TYPE")
# Convert to data frame with x/y coords
env_df <- as.data.frame(env_stack, xy = TRUE, na.rm = FALSE)

# Convert habitat to factor (and match model)
env_df$A_VEG_TYPE <- factor(env_df$A_VEG_TYPE,
                         levels = levels(unique(p_apollo_points$A_VEG_TYPE_vegetation_map)))  # must match GLM
# Step 3: Rasterize the polygon (using a specific field to assign values)

## resample with the gr_1km 
predictors <- c(gr_1km_rast,slope_rast_r,dem_rast_r,bio12_rast_r,bio15_rast_r,veg_rast_r)  # Stack the rasters
plot_rasters <- c(slope_rast_r,dem_rast_r,bio12_rast_r,bio15_rast_r,veg_rast_r)  # Stack the rasters
## plot
png("../figures/stacked_raster_with_points.png", width = 2000, height = 1500, units="px")
n <- nlyr(plot_rasters)
plot(plot_rasters, nr = 2, nc = 3)
dev.off()

## extract values
presvals <- extract(predictors, apollo_points)
# remove the ID variable
presvals <- presvals[,-1]
# setting random seed to always create the same
# random set of points for this example
set.seed(0)
backgr <- spatSample(predictors, 1000, "random", as.points=TRUE, na.rm=TRUE)
backgr$A_VEG_TYPE <- as.factor(backgr$A_VEG_TYPE)
absvals <- values(backgr)
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))

png("../figures/sdmdata_pairs.png", width = 2000, height = 1500, units="px")
pairs(sdmdata[,3:7], cex=0.8)

dev.off()

## Model fitting
m1 <- glm(pb ~ eudem_slop_3035_europe + eudem_dem_4258_europe + wc2.1_30s_bio_12 + wc2.1_30s_bio_15 + A_VEG_TYPE, data=sdmdata, family = binomial)
class(m1)
summary(m1)

#m2 = glm(pb ~ eudem_slop_3035_europe + eudem_dem_4258_europe + wc2.1_30s_bio_12 + wc2.1_30s_bio_15, data=sdmdata)
#m2

## predicts
library(predicts)
env_df$predicted_prob <- predict(m1, newdata = env_df, type = "response")

# Create a new raster layer for prediction
pred_rast <- rast(env_stack[[1]])  # use template raster
values(pred_rast) <- env_df$predicted_prob

library(RColorBrewer)
colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
png("../figures/sdm_predict_apollo.png", width = 2000, height = 1500, units="px")
plot(pred_rast, main = "Predicted Probability of Presence")
points(apollo_points)
dev.off()

### fRP
filtered_raster <- p > 0.5
masked_raster <- mask(p, filtered_raster)
cell_area <- res(p)[1] * res(p)[2]
area <- sum(!is.na(masked_raster[])) * cell_area

png("../figures/sdm_predict_apollo_0_6.png", width = 2000, height = 1500, units="px")
plot(filtered_raster,col = colors)
points(apollo_points)
dev.off()

## Model evaluation

