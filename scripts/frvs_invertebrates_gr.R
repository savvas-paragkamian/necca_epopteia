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
library(biomod2)
library(car)
library(units)
source("necca_spatial_functions.R")

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
species_samples_art17 <- read_delim("../results/species_samples_art17_private.tsv", delim="\t")

species_taxonomy <- read_delim("../results/species_gbif_taxonomy_curated.tsv",delim="\t")
iucn_art17_invert_all <- read_delim("../results/iucn_art17_invert_all.tsv", delim="\t")
iucn_art17_invert_no_tax <- iucn_art17_invert_all |>
    dplyr::select(-ends_with("Name"),submittedName)

species_dist_national_rep <- sf::st_read("../spatial_data/National report_2013_2018_shp/GR_Art17_species_distribution.shp")
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

## data preparation
apollo <- p_apollo_points |>
    st_drop_geometry() |>
    dplyr::select(decimalLatitude,decimalLongitude)


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
    geom_sf(points_convex, mapping=aes(color="gray70"), fill="NA")+
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
# add the hottest period, maybe important
# add the habitat categories which are more coarce in the categories. They are only in Natura2000
# or join some vegetation categories into one to reduce their abundance.
# or corine raster file
## add NDVI ? mean? which years? prioritize based on time availability
## raster paths


# -------------------------
# 0. Species Data P. apollo
# -------------------------
species_name <- unique(p_apollo_points$submittedName)

coords <- p_apollo_points |>
    st_drop_geometry() |> 
    dplyr::select(decimalLongitude,decimalLatitude) 

coords_rast <- vect(apollo, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

# -------------------------
# 1. Load & Resample Rasters
# -------------------------
# the 1X1 km raster that the other rasters will be resampled on.
gr_1km_terra <- file.path("../spatial_data/eea_reference_grid/gr_1km_terra.tif")
gr_1km_res <- rast(gr_1km_terra)

## individual rasters to stack
slope <- file.path("../spatial_data/EU_DEM_slope_gr/crop_eudem_slop_3035_europe.tif")
dem <- file.path("../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_4326_gr.tif")
vegetation_map <- file.path("../spatial_data/Vegetation_map_Greece/D_xabxg_VPG_60-98_GEO.tif")
corine_gr <- file.path("../spatial_data/u2018_clc2018_v2020_20u1_raster100m_gr/crop_U2018_CLC2018_V2020_20u1.tif")

individual_rasters <- c(slope,dem,vegetation_map,corine_gr)

individual_r <- lapply(individual_rasters, rast)

# the extends are different so resample is needed
lapply(individual_r, ext)

# Resample to match rows, columns, extent, and resolution
individual_res <- lapply(individual_rasters, function(f) {
                     r <- rast(f)
                     r_proj <- project(r, gr_1km_res)         # reproject if CRS differs
                     resample(r_proj, gr_1km_res, method="bilinear")  # or method="near" for categorical
})

individual_stack <- do.call(c,individual_res)

## bioclim rasters
raster_paths_all <- list.files(path = "../spatial_data",
           pattern = ".tif$",
           recursive=T,
           full.names=T)
# remove hilda and gr 

# keep only bioclimatic
bio_raster_paths <- raster_paths_all[grep("bio_",raster_paths_all,invert=F)] 

## load rasters to stack
bio_raster_stack <- rast(bio_raster_paths)

bio_raster_res <- resample(bio_raster_stack,gr_1km_res, method="bilinear")

# all rasters stack
stacked_rasters <- c(individual_stack,bio_raster_res)

# -------------------------
# 2. Crop and plot rasters
# -------------------------
# hull for crop in order to 
# reduce the area
hull_vect <- vect(points_convex)

# Crop each raster using the hull
env_stack <- crop(stacked_rasters,ext(hull_vect))

## plot
png("../figures/stacked_raster_with_points.png", width = 3000, height = 3000, units="px")
plot(env_stack)
dev.off()

# Keep only numeric layers
# Remove by name
num_stack <- env_stack[[!names(env_stack) %in% c("A_VEG_TYPE","LABEL3")]]
# Keep only categorical layers
cat_stack <- env_stack[[names(env_stack) %in% c("A_VEG_TYPE","LABEL3")]]

# -------------------------
# 2. VIF Filtering (car)
# -------------------------
# Sample random background points for VIF evaluation
set.seed(1)
bg <- spatSample(num_stack, size = 10000, method = "random",
                 na.rm = TRUE, as.points = TRUE)

# Extract values
X <- terra::extract(num_stack, bg, ID = FALSE) %>% na.omit()

# Drop non-numeric predictors
is_num <- sapply(X, is.numeric)
X <- X[, is_num, drop = FALSE]

# ------------------------
# drop highly correlated variables
# ------------------------

# function to find highly correlated variables and exclude them

find_highly_correlated <- function(cor_mat, cutoff = 0.99) {

    cor_mat[lower.tri(cor_mat, diag = TRUE)] <- 0  # keep only upper triangle
    high_corr <- which(abs(cor_mat) > cutoff, arr.ind = TRUE)
    if (nrow(high_corr) == 0) return(integer(0))  # nothing to remove
    vars_to_remove <- unique(rownames(high_corr))  # choose to drop first of each pair
    which(colnames(cor_mat) %in% vars_to_remove)
}
# correlation
cors <- cor(X, use = "pairwise.complete.obs")
cors_mat <- cors
cors[lower.tri(cors, diag = TRUE)] <- NA  # keep only upper triangle

# Convert to long format

df_cors <- as_tibble(as.data.frame(as.table(cors)))
colnames(df_cors) <- c("Var1", "Var2", "Correlation")

# plot a heatmap
cor_map_vif <- ggplot(df_cors, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1),na.value = "white", space = "Lab") +
    geom_text(aes(label = round(Correlation, 2)), size = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) +
    labs(title = "Correlation Matrix Heatmap", x = "", y = "")

ggsave("../figures/vif_correlation_matrix.png", 
       plot=cor_map_vif, 
       height = 20, 
       width = 25,
       dpi = 300, 
       units="cm",
       device="png")

# exclude the variables that have high correlation
high_cor_indices <- find_highly_correlated(cors_mat, cutoff = 0.88)

# Drop them
X_cleaned <- X[, -high_cor_indices]

# Helper function for VIF filtering
car_vif_drop <- function(X, thresh = 10) {

    X$response__ <- 1
    vars <- setdiff(names(X), "response__")
    keep <- vars
    repeat {
        mod <- lm(reformulate(keep, "response__"), data = X)
        v <- car::vif(mod)
        if (max(v) < thresh) break
        keep <- setdiff(keep, names(which.max(v)))
    }
    keep
}

selected_vars <- car_vif_drop(X_cleaned, thresh = 7)

# keep only the selected vars
preds_1km <- num_stack[[selected_vars]]

# add the categorical

preds_1km <- c(preds_1km,cat_stack)

# -------------------------
# 3. Format Data for biomod2
# -------------------------
myBiomodData <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(coords)),
  resp.xy = coords,
  resp.name = species_name,
  expl.var = stack(preds_1km),
  PA.nb.rep = 1,
  PA.nb.absences = 1000,
  PA.strategy = 'random'
)

# -------------------------
# 4. Run Single Models
# -------------------------
#myBiomodOptions <- bm_ModelingOptions("binary","default")
myBiomodOptions <- bm_ModelingOptions("binary",
                                      models = c("GLM", "RF", "ANN"),
                                      strategy="default",
                                      bm.format = myBiomodData)
# Model single models
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = 'AllModels',
                                    CV.strategy = 'random',
                                    CV.nb.rep = 2,
                                    CV.perc = 0.8,
                                    OPT.strategy = 'bigboss',
                                    var.import = 3,
                                    metric.eval = c('TSS','ROC'))
                                    # seed.val = 123)
                                    # nb.cpu = 8)
myBiomodModelOut

# Get evaluation scores & variables importance
get_evaluations(myBiomodModelOut)
var.imp <- get_variables_importance(myBiomodModelOut)

# Average across runs (3rd dimension)
var_imp_df <- var.imp |>
    group_by(algo,expl.var) |>
    summarise(Importance=mean(var.imp), .groups="keep")

colnames(var_imp_df) <- c("Model","Variable",  "Importance")

var_imp_plot <- ggplot(var_imp_df, aes(x = Model, y = Variable, fill = Importance)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Variable Importance Heatmap",
    x = "Model",
    y = "Variable",
    fill = "Importance"
  )

ggsave("../figures/biomod_var_imp_plot.png", 
       plot=var_imp_plot, 
       height = 20, 
       width = 25,
       dpi = 300, 
       units="cm",
       device="png")

# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodModelOut)
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

# Represent response curves
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[c(1:3, 12:14)],
                      fixed.var = 'min')
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = get_built_models(myBiomodModelOut)[3],
                      fixed.var = 'median',
                      do.bivariate = TRUE)

# project single model
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'Current',
                                  new.env = preds_1km,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
myBiomodProj

png("../figures/p.apollos_myBiomodProj.png", width = 3000, height = 3000, units="px")
plot(myBiomodProj)
dev.off()

proj_rast <- get_predictions(myBiomodProj)

## Response curves
## ROC curves
# -------------------------
# 5. Run Ensemble Models
# -------------------------
# Model ensemble models
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      em.algo = c('EMmean', 'EMcv', 'EMci', 'EMmedian', 'EMca', 'EMwmean'),
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.7),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 3,
                                      EMci.alpha = 0.05,
                                      EMwmean.decay = 'proportional')
myBiomodEM

# Get evaluation scores & variables importance
get_evaluations(myBiomodEM)
get_variables_importance(myBiomodEM)

# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = myBiomodEM, group.by = 'full.name')
bm_PlotEvalBoxplot(bm.out = myBiomodEM, group.by = c('full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'full.name', 'full.name'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'merged.by.run'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'expl.var', 'merged.by.run'))

# Represent response curves
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[c(1, 6, 7)],
                      fixed.var = 'min')
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM)[7],
                      fixed.var = 'median',
                      do.bivariate = TRUE)

# Project ensemble models (from single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                             bm.proj = myBiomodProj,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')


myBiomodProj
plot(myBiomodProj)

# -------------------------
# 6. Ensemble Modeling
# -------------------------

# Project ensemble models (from single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                             bm.proj = myBiomodProj,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
                                             
# Project ensemble models (building single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                             proj.name = 'CurrentEM',
                                             new.env = myExpl,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
myBiomodEMProj
plot(myBiomodEMProj)



myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = "all",
  em.by = "all",
  eval.metric = c("TSS", "ROC"),
  eval.metric.quality.threshold = c(0.7, 0.8)
)

# -------------------------
# 7. Projection on 1x1 km Grid
# -------------------------
myBiomodProjection <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = stack(preds_1km),
  proj.name = "current_1km",
  selected.models = "all",
  binary.meth = "TSS",
  build.clamping.mask = FALSE
)

# -------------------------
# 8. Save Output
# -------------------------
proj_rast <- get_predictions(myBiomodProjection)  # RasterStack
writeRaster(proj_rast, "SDM_current_1km.tif", overwrite = TRUE)





# 3) Sample background points in M
set.seed(1)
preds <- env_stack
bg <- spatSample(preds,
                 size = 20000, method = "regular", na.rm = TRUE, as.points = TRUE, warn = FALSE)
bg <- bg[hull_vect]  # keep only inside M
## plot
png("../figures/stacked_raster_sampled.png", width = 3000, height = 3000, units="px")
plot(bg)
dev.off()

# 4) Extract env values for VIF

X <- terra::extract(preds, bg, ID = FALSE)
X <- na.omit(as.data.frame(X))

# 5) Run VIF
v <- vifstep(X, th = 10)   # or th = 5
v@results
vars_kept <- v@results$Variables

# 6) Keep only selected predictors
preds_sel <- preds[[vars_kept]]


# Convert to data frame with x/y coords
env_df <- as.data.frame(env_stack, xy = TRUE, na.rm = FALSE)

# Convert habitat to factor (and match model)
env_df$A_VEG_TYPE <- factor(env_df$A_VEG_TYPE,
                         levels = levels(unique(p_apollo_points$A_VEG_TYPE_vegetation_map)))  # must match GLM

env_df_species <- env_df |> 
    filter(!is.na(x))

#################### VIF analysis 
###

# Step 3: Rasterize the polygon (using a specific field to assign values)

## resample with the gr_1km 
#predictors <- c(gr_1km_rast,slope_rast_r,dem_rast_r,bio12_rast_r,bio15_rast_r,veg_rast_r)  # Stack the rasters
predictors <- env_stack

## extract values
presvals <- extract(env_stack, apollo_points)
# setting random seed to always create the same
# random set of points for this example
set.seed(0)

# minimum convex hull or bbox to use as a background in order to exclude non present areas
# 
backgr <- spatSample(predictors, 1000, "random", as.points=TRUE, na.rm=TRUE)
backgr$A_VEG_TYPE <- as.factor(backgr$A_VEG_TYPE)
absvals <- values(backgr)
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))

png("../figures/sdmdata_pairs.png", width = 2000, height = 1500, units="px")
pairs(sdmdata[,3:7], cex=0.8)

dev.off()

### 75% train and 25% predict to evaluate the performance

## Model fitting
m1 <- glm(pb ~ eudem_slop_3035_europe + eudem_dem_4258_europe + wc2.1_30s_bio_12 + wc2.1_30s_bio_15 + A_VEG_TYPE, data=sdmdata, family = binomial)
class(m1)
summary(m1)

### focus on maxend.
### Random forest? 
### future ensemble forcasting

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
p <- pred_rast
filtered_raster <- p > 0.5
masked_raster <- mask(p, filtered_raster)
cell_area <- res(p)[1] * res(p)[2]
area <- sum(!is.na(masked_raster[])) * cell_area

r_template <- rast(locations_10_grid_samples, resolution = 100, crs = st_crs(locations_10_grid_samples)$wkt)

r <- rasterize(vect(locations_10_grid_samples),r_template ,field = "n_samples")

png("../figures/sdm_predict_apollo_0_6.png", width = 2000, height = 1500, units="px")
plot(r,col = "yellow")
plot(filtered_raster,col = colors)
points(apollo_points)
dev.off()

## Model evaluation

