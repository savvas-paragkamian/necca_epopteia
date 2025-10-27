#!/usr/bin/env Rscript

## Script name: frvs_invertebrates_gr.R
##
## Purpose of script: Calculate the Favourite Reference Values
## of the invertebrates in Greece
##
## Author: Savvas Paragkamian, Christina Kassara
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

greece_regions <- sf::st_read("../spatial_data/gadm41_GRC_shp/gadm41_GRC_2.shp") |>
    st_transform(crs="EPSG:3035")

gr_1km <- sf::st_read("../spatial_data/eea_reference_grid/gr_1km.shp") 

gr_10km <- sf::st_read("../spatial_data/eea_reference_grid/gr_10km.shp") 

N2000_v32 <- sf::st_read("../spatial_data/N2000_spatial_GR_2021_12_09_v32/N2000_spatial_GR_2021_12_09_v32.shp")

N2000_v32_wgs <- st_transform(N2000_v32,4326)

########################### Load Species Data ###########################
### Species occurrences enriched ######
species_samples_art17 <- read_delim("../results/species_samples_presence_final_private.tsv", delim="\t")

            
species_taxonomy <- read_delim("../results/species_gbif_taxonomy_curated.tsv",delim="\t")
iucn_art17_invert_all <- read_delim("../results/iucn_art17_invert_all.tsv", delim="\t")
iucn_art17_invert_no_tax <- iucn_art17_invert_all |>
    dplyr::select(-ends_with("Name"),submittedName)

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

parnassius_dist <- sf::st_read("../data/Parnassius apollo AP 2019/AP_Papollo_Distribution_LAEA.shp") 

p_apollo_points_wgs <- species_art17_spatial |>
    filter(species=="Parnassius apollo") |>
    filter(includeDistribution==TRUE | datasetName=="GBIF") |>
    filter(datasetName!="Action Plan 2019")
    #filter(!c(datasetName %in% c("Action Plan 2019", "DistrMap_2013_2018")))

p_apollo_points <- p_apollo_points_wgs |>
    st_transform(3035)

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
                     "Action Plan 2019"="#B31319",
                     "E1X_MDPP_2014_2024"="#FDF79C",
                     "DistrMap_2013_2018"="black",
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
    geom_sf(N2000_v32, mapping=aes(fill=SITETYPE),
            alpha=0.8,
            #colour="transparent",
            na.rm = F,
            show.legend=T) +
#    geom_sf(points_convex, mapping=aes(color="gray70"), fill="NA")+
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

# --------------------------------------------------------- #
###################### Modeling FRV ########################
# --------------------------------------------------------- #

# -------------------------
# 0. Species Data P. apollo
# -------------------------
species_name <- unique(p_apollo_points$submittedName)

## data preparation
points_convex <- p_apollo_points |>
    st_union() |>
    st_convex_hull() |>
    st_buffer(dist = 1000)

st_write(points_convex,"../results/p_apollo/geospatial/convex_hull.shp",overwrite=T,append=FALSE)

coords <- st_coordinates(p_apollo_points)

# ---------
# convex hull
# __________
# hull for crop in order to 
# reduce the area
hull_vect <- vect(points_convex)

# ------------------------------------
# 1. Load & Resample Numerical Rasters
# ------------------------------------

## individual rasters to stack
slope_f <- file.path("../spatial_data/EU_DEM_slope_gr/crop_eudem_slop_3035_europe.tif")
dem_f <- file.path("../spatial_data/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_3035_europe.tif")

# slope load
slope <- rast(slope_f)
# dem load
dem <- rast(dem_f)

individual_rasters <- c(slope,dem)

## bioclim rasters
raster_paths_all <- list.files(path = "../spatial_data",
           pattern = ".tif$",
           recursive=T,
           full.names=T)
# remove hilda and gr 

# keep only bioclimatic
bio_raster_paths <- raster_paths_all[grep("bio_",raster_paths_all,invert=F)] 

## load rasters to stack which are in WGS84
bio_raster_stack_w <- rast(bio_raster_paths)

# change the projection to "EPSG:3035"
bio_raster_stack <- project(bio_raster_stack_w, "EPSG:3035")

#bio_raster_stack <- project(bio_raster_stack_w, gr_1km_res)

#bio_raster_res <- resample(bio_raster_stack,gr_1km_res, method="bilinear")

# ---------------------------------------------
# 2. Load Classify Categorical Rasters
# ---------------------------------------------


#gr_1km_res <- rast("../spatial_data/eea_reference_grid/EEA_grid_id_1km.tif")

# corine land cover in raster format
corine_gr <- rast("../spatial_data/u2018_clc2018_v2020_20u1_raster100m_gr/crop_U2018_CLC2018_V2020_20u1.tif")

# -------------------------------------------
# first reclassify based on manual curation
# according to relavant categories resolution
# -------------------------------------------
# clc
corine <- as.factor(corine_gr)

# vegetation map
#vegetation_map_r <- rast(vegetation_map)

# summary of the categorical variables
freq_corine <- freq(corine_gr)

#freq_vegtype <- freq(vegetation_map_r)

# extract the points
# corine
corine_points <- terra::extract(corine_gr,p_apollo_points)

# vegetation
#veg_points <- terra::extract(vegetation_map_r,p_apollo_points) 

cat_points <- corine_points |> 
#    left_join(veg_points) |>
    count(LABEL3)

#points_summary_veg <- veg_points |>
#    group_by(A_VEG_TYPE) |>
#    summarise(n_points=n())

# label 3 corine
points_summary_label3 <- corine_points |>
    group_by(LABEL3) |>
    summarise(n_points=n())

# after this analysis we move to re-classification
# of the corine raster in order to reduce categories
# and keep the relavant ones for P. apollo

# manual matrix the maps to the new categories
rcl <- matrix(c(
                1,1,
                2,1,
                3,1,
                4,1,
                5,1,
                6,1,
                7,1,
                8,1,
                9,1,
               10,1,
               11,1,
               12,2,
               13,2,
               14,2,
               15,2,
               16,2,
               17,2,
               18,2,
               19,2,
               20,2,
               21,2,
               22,2,
               23,23,
               24,24,
               25,25,
               26,26,
               27,27,
               28,28,
               29,29,
               30,30,
               31,31,
               32,32,
               33,33,
               34,5,
               35,5,
               36,5,
               37,5,
               38,5,
               39,5,
               40,5,
               41,5,
               42,5,
               43,5,
               44,NA
), ncol=2, byrow=TRUE)

# now reclassify
corine_r <- classify(corine,rcl)

#values(corine_r)[values() < 0] = NA

# summary of the re classified corine
corine_r_freq <- freq(corine_r)

# ---------------------------------------------
# 3. Crop and resample Rasters to 1X1 km Res
# ---------------------------------------------


bio1 <- crop(bio_raster_stack[[1]],ext(dem))

png("../figures/rasters_gr.png",
    width = 6000,
    height = 1500,
    res=300,
    units="px")
par(mfrow = c(1,4),        # 1 row, 4 columns
    oma = c(1,1,1,1),      # outer margins
    mar = c(1,1,1,1),      # inner margins
    cex = 1.2)             # text size
plot(corine, legend = FALSE, main = "CORINE Land Cover")
plot(dem, legend = FALSE, main = "Digital Elevation Model")
plot(slope, legend = FALSE, main = "Slope")
plot(bio1, legend = FALSE, main = "Bioclim Layer 1")
dev.off()

# -----------------
# bounding box 
# using the template of gr 1km reference grid
# to create a raster to use in the models
# -----------------

# load the grid 1km 
shp <- gr_1km

# manual bbox of the area of interest
xmin = 5200000
xmax = 5570000
ymin = 1720000
ymax = 2180000

bb <- st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = st_crs(shp))

# find the intersecting polygons
sel <- st_intersects(shp, st_as_sfc(bb), sparse = FALSE)[,1]

# geometries remain intact
shp_bbox <- shp[sel, ]

# create a polygon to use as mask with an extent
# and havind the overlap with the 1km EEA reference grid
e <- terra::ext(xmin, xmax, ymin, ymax)

p <- terra::as.polygons(e, crs="EPSG:3035")

#bbox_0 <- terra::crop(p, terra::ext(shp))
bbox <- terra::crop(p, e)

# the 1X1 km raster that the other rasters will be resampled on.
bbox_grid <- terra::rasterize(bbox,terra::rast(e,resolution = 1000)) # 1km res

# assign the CRS manually
crs(bbox_grid) <- "EPSG:3035"

# export to validate in gis GUI apps
writeRaster(bbox_grid, "../results/p_apollo/geospatial/bbox_grid.tif", overwrite = TRUE)

# -----------------
## Crop and resample
## each raster using the coun.rast
# -----------------

bio_c <- crop(bio_raster_stack, bbox_grid)

bio_raster_stack_r <- resample(bio_raster_stack, bbox_grid, method="average", threads=TRUE)

dem_c <- crop(dem, bbox_grid)

dem_r <- resample(dem_c, bbox_grid, method="average", threads=TRUE)

writeRaster(dem_c, "../results/p_apollo/geospatial/dem_c.tif", overwrite = TRUE)
writeRaster(dem_r, "../results/p_apollo/geospatial/dem_r.tif", overwrite = TRUE)

slope_c <- crop(slope, bbox_grid)

slope_r <- resample(slope_c, bbox_grid, method="average", threads=TRUE)
writeRaster(slope_r, "../results/p_apollo/geospatial/slope_r.tif", overwrite = TRUE)


# Corine
## crop doesn\'t work with factor rasters. 
## so i use the coun.crop
corine_rc <- crop(corine_r,bbox_grid)

writeRaster(corine_rc, "../results/p_apollo/geospatial/corine_rc_crop.tif", overwrite = TRUE)

#corine_rcr <- resample(corine_rc, bbox_grid, method="mode", threads=TRUE)

#corine_rcr_m <- resample(corine_rc, bbox_grid, method="max", threads=TRUE)

#writeRaster(corine_rcr, "../spatial_data/corine_rct_mode.tif", overwrite = TRUE)

#writeRaster(corine_rcr_m, "../spatial_data/corine_rct_max.tif", overwrite = TRUE)

# --------------------
# analyse the overlap of corine categories
# per 1km cell
# -------------------

library(exactextractr)

# assign values to raster 
x <- setValues(bbox_grid, 1:ncell(bbox_grid))

# make vector
bbox_grid_s <- terra::as.polygons(x)

bbox_grid_sf <- sf::st_as_sf(bbox_grid_s)

st_write(bbox_grid_sf, "../results/p_apollo/geospatial/bbox_grid_sf.shp",append=TRUE)

# extract the grid ids that the p apollo has points
p_apollo_bbox <- terra::extract(x,p_apollo_points) 

p_apollo_corine <- terra::extract(corine_rc, p_apollo_points)

p_apollo_spatial <- p_apollo_bbox |> left_join(p_apollo_corine) |>
    mutate(layer=as.character(layer)) |>
    group_by(layer) |>
    summarise(p_apollo=n(),LABEL3=str_c(LABEL3,collapse=","))

colnames(p_apollo_spatial) <- c("rowname","p_apollo","LABEL3")

#write_delim(extract_count, "../results/extract_count.tsv", delim="\t")

# ---------------------
# tabulate area of bbox grid with the re classified corine
# --------------------

# exact_extract is a very fast function from the exactextractr package
extract_frac = exact_extract(corine_rc,bbox_grid_sf,fun="frac")

extract_frac_t <- as_tibble(extract_frac) |>
    rownames_to_column()

write_delim(extract_frac_t, "../results/p_apollo/geospatial/extract_frac.tsv", delim="\t")

# make long format to find the dominant classes

extract_frac_l <- extract_frac_t |>
    pivot_longer(-rowname,names_to="categories", values_to="frac")

# greater than 0.5 is dominant, and keep the value of frac,
# if it is not greater, then check if sum == 0, meaning
# this cell has no value in the corine categories, and assign the value 200.
# In all other cases, set -1
extract_frac_dom <- extract_frac_l |>
    group_by(rowname) |>
    mutate(sum=sum(frac),
           max=max(frac),
           max_variable = categories[which.max(frac)]) |>
    mutate(dominant = if_else( sum==0,200, if_else(max>0.5,max,-1))) |>
    #mutate(cat = if_else( sum==0,"zero", if_else(max>0.5,frac,"mixed"))) |>
    ungroup() |>
    mutate(class=if_else(dominant==200,
                         "zero",
                         if_else(dominant==-1,"mixed","dominant"))) 

# summarise based on the cell
extract_frac_dom_dist <- extract_frac_dom |>
    distinct(rowname,sum,dominant,class,max_variable)

# make a summary of the cells according to their category
extract_frac_dom_sum <- extract_frac_dom |>
    distinct(rowname,dominant,class) |> 
    group_by(class) |>
    summarise(n=n())

# join with the original table
# manually create new columns based on manually curated 
# corine categories
# clc (23,24,25,26,27,29,31,32)
#Create 3 new rasters: Sum percentage of clc (23,24,25,26,27,29,31,32) / Percentage of clc 26 / Sum percentage of clc (23,24,25,27,29,31,32). Caution: if sum==0, then value=NA in new rasters

extract_frac_dom_s <- extract_frac_t |>
    left_join(extract_frac_dom_dist) |>
    left_join(p_apollo_spatial) |>
    ungroup() |>
    mutate(PercSuitClassAll=if_else(sum!=0,frac_23+frac_24+frac_25+frac_26+frac_27+frac_29+frac_31+frac_32,NA)) |>
    mutate(PercSuitClassElse=if_else(sum!=0,frac_23+frac_24+frac_25+frac_27+frac_29+frac_31+frac_32,NA)) |>
    mutate(PercFrac_26=if_else(sum!=0,frac_26,NA))

write_delim(extract_frac_dom_s,"../results/p_apollo/geospatial/extract_frac_dom_s.tsv",delim="\t")

# -------------------------
# create the rasters based 
# on the vectors that were 
# created before
# -------------------------

#### PercSuitClassAll_m
PercSuitClassAll_m <- as.matrix(data.frame(as.numeric(extract_frac_dom_s$rowname),as.numeric(extract_frac_dom_s$PercSuitClassAll)))

PercSuitClassAll_r <- classify(x,PercSuitClassAll_m)
names(PercSuitClassAll_r) <- "PercSuitClassAll"

writeRaster(PercSuitClassAll_r, "../results/p_apollo/geospatial/PercSuitClassAll_r.tif", overwrite = TRUE)


#### PercSuitClassElse
PercSuitClassElse_m <- as.matrix(data.frame(as.numeric(extract_frac_dom_s$rowname),as.numeric(extract_frac_dom_s$PercSuitClassElse)))

PercSuitClassElse_r <- classify(x,PercSuitClassElse_m)
names(PercSuitClassElse_r) <- "PercSuitClassElse"

writeRaster(PercSuitClassElse_r, "../results/p_apollo/geospatial/PercSuitClassElse_r.tif", overwrite = TRUE)

#### PercFrac_26
PercFrac_26_m <- as.matrix(data.frame(as.numeric(extract_frac_dom_s$rowname),as.numeric(extract_frac_dom_s$PercFrac_26)))

PercFrac_26_r <- classify(x,PercFrac_26_m)
names(PercFrac_26_r) <- "PercFrac_26"

writeRaster(PercFrac_26_r, "../results/p_apollo/geospatial/PercFrac_26_r.tif", overwrite = TRUE)



# numeric stack
stacked_rasters_n <- c(PercFrac_26_r,PercSuitClassElse_r,PercSuitClassAll_r,bio_raster_stack_r, dem_r, slope_r)

stacked_rasters_examp <- c("PercSuitClassAll",
                        "eudem_dem_3035_europe",
                        "eudem_slop_3035_europe",
                        "wc2.1_30s_bio_1")

rasters_gr <- stacked_rasters_n[[stacked_rasters_examp]]

# plot example
png("../figures/stacked_raster_example.png",
    width = 4000,
    height = 1500,
    res=300,
    units="px")
#par(cex = 5,    # global expansion for text
#    oma = c(4,4,4,4),   # outer margins
#    mar = c(5,5,4,2))   # inner margins
plot(rasters_gr,
     nr=1,
     nc=4
)
dev.off()

# ----------------------------------
# save rasters to ascii
# ----------------------------------
names_stacked <- names(stacked_rasters_n)

for (i in seq_along(names_stacked)){
    
    name <- names_stacked[i]

    writeRaster(
      stacked_rasters_n[[i]],
      filename = paste0("../results/p_apollo/geospatial/ascii/num_stack_",name,".asc"),
      #filetype = "ascii",      # tells terra to write ESRI ASCII Grid
      NAflag   = -9999,        # NA value written in the ASCII file
      overwrite = TRUE
    )
}



# ----------------------------------
# Crop each raster using the hull
# ----------------------------------
num_stack <- crop(stacked_rasters_n,ext(hull_vect))


# extract the values
p_apollo_stack_a <- terra::extract(stacked_rasters_n,p_apollo_points,xy=TRUE)

names(p_apollo_stack_a) <- sub("wc2.1_30s_","",names(p_apollo_stack_a))

names(p_apollo_stack_a) <- sub("eudem_slop_3035_europe","slop",names(p_apollo_stack_a))
names(p_apollo_stack_a) <- sub("eudem_dem_3035_europe","dem",names(p_apollo_stack_a))


p_apollo_stack <- terra::extract(num_stack,p_apollo_points)


write_delim(p_apollo_points,"../results/p_apollo/p_apollo_points.tsv",delim="\t")

write_delim(p_apollo_stack_a,"../results/p_apollo/p_apollo_stack.csv",delim=",")


p_apollo_stack_a_shp <- p_apollo_stack_a |>
    st_as_sf(coords=c("x","y"),
             remove=F,
             crs="EPSG:3035") 


st_write(p_apollo_stack_a_shp, "../results/p_apollo/p_apollo_stack_a_shp.shp",delete_dsn = TRUE)

# Keep only numeric layers
# Remove by name
#num_stack <- env_stack[[!names(env_stack) %in% c("A_VEG_TYPE","LABEL3")]]
# Keep only categorical layers
#cat_stack <- env_stack[[names(env_stack) %in% c("A_VEG_TYPE","LABEL3")]]

# -------------------------
# 4. VIF Filtering (car)
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

ggsave("../figures/correlation_matrix.png", 
       plot=cor_map_vif, 
       height = 20, 
       width = 25,
       dpi = 300, 
       units="cm",
       device="png")

# exclude the variables that have high correlation
high_cor_indices <- find_highly_correlated(cors_mat, cutoff = 0.7)

keep <- c(
          "eudem_dem_3035_europe",
          "eudem_slop_3035_europe",
          "wc2.1_30s_bio_3",
          "wc2.1_30s_bio_4",
          "wc2.1_30s_bio_8",
          "wc2.1_30s_bio_9",
          "wc2.1_30s_bio_13",
          "wc2.1_30s_bio_14",
          #"PercFrac_26",
          #"PercSuitClassElse",
          "PercSuitClassAll"
)

# Drop them
X_cleaned <- X[, keep]

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
    list(v,keep)
}

selected_vars <- car_vif_drop(X_cleaned, thresh = 5) # keep 5 which is commonly used.

# keep only the selected vars
preds_1km_bbox <- num_stack[[selected_vars[[2]]]]

# ------------------------
# Plot
# ------------------------

## plot 
png("../figures/stacked_raster_with_points.png",
    width = 4000,
    height = 2500,
    res=300,
    units="px")
#par(cex = 5,    # global expansion for text
#    oma = c(4,4,4,4),   # outer margins
#    mar = c(5,5,4,2))   # inner margins
plot(preds_1km_bbox,
     nr=2,
     nc=4
)
dev.off()

# mask (crop with vector) with convex hull
preds_1km <- mask(preds_1km_bbox, hull_vect)

# -------------------------
## check for NA before running the model
# -------------------------
which(is.na(terra::extract(preds_1km, p_apollo_points)))

# -------------------------
# function for models
# -------------------------

models_function <- function(myBiomodData,model_name, preds) {
    
    model_name <- as.character(model_name)
    print(model_name)
    
    preds_1km <- preds
    
    # ---------------------------------
    # plot the rasters
    # ---------------------------------
    png(paste0("../figures/",model_name,"_stacked_raster_with_points.png"),
        width = 4000,
        height = 2500,
        res=300,
        units="px")
    #par(cex = 5,    # global expansion for text
    #    oma = c(4,4,4,4),   # outer margins
    #    mar = c(5,5,4,2))   # inner margins
    plot(preds_1km,
         nr=2,
         nc=4
    )
    dev.off()
    
    #stop("stop here")
    #}
    # -------------------------
    # Cross validation
    # -------------------------
    
    #cv.r <- bm_CrossValidation(bm.format = myBiomodData,
    #                           strategy = "random",
    #                           nb.rep = 3,
    #                           k = 0.8)
    #
    ## k-fold selection
    #cv.k <- bm_CrossValidation(bm.format = myBiomodData,
    #                           strategy = "kfold",
    #                           nb.rep = 2,
    #                           k = 3)
    #
    #
    # -------------------------
    # 6. Run Single Models
    # -------------------------
    ## Steven J. Phillips, Miroslav Dudík, Robert E. Schapire. [Internet] Maxent software for modeling species niches and distributions (Version 3.4.1). Available from url: http://biodiversityinformatics.amnh.org/open_source/maxent/. Accessed on 2025-7-31.
    #myBiomodOptions <- bm_ModelingOptions("binary","default")
    
    # GLM could not finish
    #my_models <- c("GLM", "RF", "MAXENT")
    
    #my_models <- c("MAXENT","MAXNET")
    #my_models <- c("MAXNET")
    my_models <- c("RF")
    #my_models <- c("GLM", "RF")

    #user.MAXENT_official <- list('for_all_datasets' = list(
    #                                                       betamultiplier = 0.2,
    #                                                       linear = TRUE,
    #                                                       quadratic = TRUE,
    #                                                       jackknife = TRUE,
    #                                                       product = TRUE,
    #                                                       lq2lqptthreshold = 0,
    #                                                       outputformat="cloglog",
    #                                                       hingethreshold = 0,
    #                                                       threshold = FALSE,
    #                                                       hinge = TRUE
    #                                                       ))

    #opts <- bm_ModelingOptions(RF = list( ntree = 1500, mtry = floor(sqrt(ncol(preds))),nodesize = 5, sampsize = NULL, replace = TRUE ))#user_val <- list(MAXENT.binary.MAXENT.MAXENT=user.MAXENT_official)
    
    #myBiomodOptions <- bm_ModelingOptions("binary",
    #                                      models = my_models,
    #                                      strategy="default",
    ##                                      user.val = user.MAXENT_official,
    #                                      bm.format = myBiomodData)

    # Model single models
    
    myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                        modeling.id = "P_apollo",
                                        models = my_models,
                                        CV.strategy = 'kfold',
                                        CV.k = 5,
                                        CV.nb.rep = 10,
                                        CV.perc = 0.7, # data percent for calibration
                                        OPT.strategy = 'bigboss',
    #                                    OPT.user = myBiomodOptions,
                                        var.import = 1,
                                        #nb.cpu = 4,
                                        do.progress=T,
                                        metric.eval = c("ROC","TSS"))
                                        # seed.val = 123)
    
    PlotEvalBoxplot <- bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
    
    
    #PlotResponseCurves <- bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
    #                      models.chosen = grep("allRun",get_built_models(myBiomodModelOut), value=T),
    #                      fixed.var = 'median',
    #                      do.bivariate = TRUE)
    
    #ggsave(paste0("../figures/",model_name,"_PlotResponseCurves.png"), 
    #       plot=PlotResponseCurves$plot, 
    #       height = 100, 
    #       width = 100,
    #       dpi = 300, 
    #       units="cm",
    #       device="png")
    
    # Get evaluation scores & variables importance
    var.eval <- get_evaluations(myBiomodModelOut)
    write_delim(var.eval,paste0("../results/",model_name,"_var.eval.tsv"),delim="\t")
    
    var.imp <- get_variables_importance(myBiomodModelOut)
    write_delim(var.imp,paste0("../results/",model_name,"_var.imp.tsv"),delim="\t")
    
    #---------------------------------------------------
    # Plots 
    # Represent evaluation scores & variables importance
    # --------------------------------------------------
    
    var.eval_g <- ggplot() +
        geom_boxplot(data=var.eval,
                     mapping=aes(x=metric.eval,y=validation)) +
      theme_bw() +
      facet_wrap(~algo)
    
    
    ggsave(paste0("../figures/",model_name,"_biomod_var_eval_plot.png"), 
           plot=var.eval_g, 
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")
    
    # Average across runs (3rd dimension)
    var_imp_df <- var.imp |>
        group_by(algo,expl.var) |>
        summarise(Importance=mean(var.imp), .groups="keep") |>
        mutate(Importance=round(Importance,digits=2))
    
    colnames(var_imp_df) <- c("Model","Variable",  "Importance")
    
    var_imp_plot <- ggplot() +
      geom_tile(var_imp_df,mapping=aes(x = Model, y = Variable, fill = Importance), color = "gray90") +
      geom_text(var_imp_df,mapping=aes(x = Model, y = Variable, label = Importance),
                size=4) +
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
    
    ggsave(paste0("../figures/",model_name,"_biomod_var_imp_plot.png"), 
           plot=var_imp_plot, 
           height = 20, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")
    
    # --------------------------------------------------
    # Represent response curves
    # --------------------------------------------------
    
    #response_c <- bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
    #                      models.chosen = get_built_models(myBiomodModelOut),
    #                      new.env = get_formal_data(myBiomodModelOut, "expl.var"),
    #                      show.variables = get_formal_data(myBiomodModelOut, "expl.var.names"),
    #                      fixed.var = 'median',
    #                      do.bivariate = F)
    #
    # --------------------------------------------------
    # Projections
    # --------------------------------------------------
    
    # project single model
    myBiomodProj <- BIOMOD_Projection(
      bm.mod   = myBiomodModelOut,           # result from BIOMOD_Modeling()
      new.env           = preds_1km,
      proj.name         = "prob_current",
      selected.models   = get_built_models(myBiomodModelOut),  # or a subset like grep("RF", ...)
      binary.meth       = NULL,                # we want continuous probs, not binary
      compress          = "xz",
      output.format     = ".img",              # or ".img" ; for GeoTIFF see below
      build.clamping.mask = TRUE,              # optional: flags extrapolation
      do.stack          = TRUE,
      keep.in.memory    = TRUE
    )
    
    # --------------------------------------------------
    # Ensemble
    # --------------------------------------------------
    
    myBiomodEM <- BIOMOD_EnsembleModeling(
                                        bm.mod = myBiomodModelOut,
                                        models.chosen = get_built_models(myBiomodModelOut),
                                        em.by           = "algo",             # <- key: one model per algorithm
                                        em.algo = c('EMmean'),
                                        metric.select = c('ROC'),
                                        #metric.select.thresh = c(0.7),
                                        var.import = 3
                                        )
    
    # Project ensemble models (from single projections)
    myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                                 bm.proj = myBiomodProj,
                                                 models.chosen = 'all',
                                                 metric.binary = 'all',
                                                 metric.filter = 'all')
    
    ## get prediction
    proj_rast <- get_predictions(myBiomodEMProj)
    writeRaster(proj_rast, paste0("../results/p_apollo/geospatial/",model_name,"_SDM_prediction_current_1km.tif"), overwrite = TRUE)
    
    ### plot
    png(paste0("../figures/",model_name,"_p.apollos_myBiomodProj.png"),
        width = 4000,
        height = 2500,
        res=300,
        units="px")
    plot(proj_rast)
    dev.off()
    
    # response ensemble
    
    response_e <- bm_PlotResponseCurves(bm.out = myBiomodEM, 
                          models.chosen = get_built_models(myBiomodEM),
                          new.env = get_formal_data(myBiomodEM, "expl.var"),
                          show.variables = get_formal_data(myBiomodEM, "expl.var.names"),
                          fixed.var = 'median',
                          do.bivariate = F)
    #i <- 0
    #for (p in plots) {
    #    if (!inherits(p, "ggplot")) next
    #    i <- i + 1
    #    ggplot2::ggsave(sprintf(paste0("../figures/",model_name,"_rcurves/rcurve_%03d.pdf"), i),
    #              plot = p,
    #              device = grDevices::cairo_pdf, width = 10, height = 7, limitsize = FALSE)
    #}
    ggsave(paste0("../figures/",model_name,"_biomod_response_ensemble.png"), 
           plot=response_e$plot, 
           height = 30, 
           width = 25,
           dpi = 300, 
           units="cm",
           device="png")

    
    # return the model
    myBiomodModelOut

############### END OF FUNCTION ###############
}
############### END OF FUNCTION ###############

#
#################################################################

models_out <- list()

#################################################################
#### model 1
#################################################################

keep <- c(
          "eudem_dem_3035_europe",
          "eudem_slop_3035_europe",
          "wc2.1_30s_bio_3",
          "wc2.1_30s_bio_4",
          "wc2.1_30s_bio_8",
          "wc2.1_30s_bio_9",
          "wc2.1_30s_bio_13",
          "wc2.1_30s_bio_14",
          #"PercFrac_26",
          #"PercSuitClassElse",
          "PercSuitClassAll"
)

preds_model1 <- num_stack[[keep]]

preds <- stacked_rasters_n[[keep]]

# prefix
model_name <- "model1"

# disk method
PA.d <- bm_PseudoAbsences(resp.var = rep(1, nrow(coords)),
                          expl.var = preds_model1,
                          nb.absences = 10000,
                          strategy = 'random',
                          dist.min = 1000)


# data
myBiomodData <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(coords)),
  resp.xy = coords,
  resp.name = species_name,
  expl.var = preds_model1,
  PA.dist.min = 1000,
  PA.nb.rep = 1,
  PA.nb.absences = 10000,
  seed.val=123,
  PA.strategy = 'random'
)

models_out[[model_name]] <- models_function(myBiomodData,model_name,preds_model1)


xy_all <- get_formal_data(myBiomodData, "coord")


#################################################################
#### model 2 
#################################################################

#Model 2 (bio2 + all clc)
#	Slope
#	Dem
#	Bio2
#	Bio4
#	Bio8
#	Bio9
#	Bio13
#	Bio14
#	Sum percentage of clc (23,24,25,26,27,29,31,32)

keep <- c(
          "eudem_dem_3035_europe",
          "eudem_slop_3035_europe",
          "wc2.1_30s_bio_2",
          "wc2.1_30s_bio_4",
          "wc2.1_30s_bio_8",
          "wc2.1_30s_bio_9",
          "wc2.1_30s_bio_13",
          "wc2.1_30s_bio_14",
          #"PercFrac_26",
          #"PercSuitClassElse",
          "PercSuitClassAll"
)

preds_model2 <- num_stack[[keep]]
preds <- stacked_rasters_n[[keep]]

# prefix
model_name <- "model2"

#preds_model2 <- mask(preds_model2, hull_vect)

myBiomodData <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(coords)),
  resp.xy = coords,
  resp.name = species_name,
  expl.var = preds_model2,
  PA.nb.rep = 1,
  PA.nb.absences = 10000,
  PA.strategy = 'random'
)

models_out[[model_name]] <- models_function(myBiomodData,"model2",preds)

#################################################################
#### model 3
#################################################################

#Model 3 (bio3 + clc26)
#1.	Slope
#2.	Dem
#3.	Bio3
#4.	Bio4
#5.	Bio8
#6.	Bio9
#7.	Bio13
#8.	Bio14
#9.	Percentage of clc 26

keep <- c(
          "eudem_dem_3035_europe",
          "eudem_slop_3035_europe",
          "wc2.1_30s_bio_3",
          "wc2.1_30s_bio_4",
          "wc2.1_30s_bio_8",
          "wc2.1_30s_bio_9",
          "wc2.1_30s_bio_13",
          "wc2.1_30s_bio_14",
          "PercFrac_26"
          #"PercSuitClassElse",
          #"PercSuitClassAll"
)

# keep only the selected rasters
preds_model3 <- num_stack[[keep]]
preds <- stacked_rasters_n[[keep]]

# prefix
model_name <- "model3"

#preds_model3 <- mask(preds_model3, hull_vect)

myBiomodData <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(coords)),
  resp.xy = coords,
  resp.name = species_name,
  expl.var = preds_model3,
  PA.nb.rep = 1,
  PA.nb.absences = 10000,
  PA.strategy = 'random'
)

models_out[[model_name]] <- models_function(myBiomodData,model_name, preds_model3)

#################################################################
#### model 4
#################################################################

#Model 4 (bio2 + clc26)
#1.	Slope
#2.	Dem
#3.	Bio2
#4.	Bio4
#5.	Bio8
#6.	Bio9
#7.	Bio13
#8.	Bio14
#9.	Percentage of clc 26

keep <- c(
          "eudem_dem_3035_europe",
          "eudem_slop_3035_europe",
          "wc2.1_30s_bio_2",
          "wc2.1_30s_bio_4",
          "wc2.1_30s_bio_8",
          "wc2.1_30s_bio_9",
          "wc2.1_30s_bio_13",
          "wc2.1_30s_bio_14",
          "PercFrac_26"
          #"PercSuitClassElse",
          #"PercSuitClassAll"
)

preds_model4 <- num_stack[[keep]]
preds <- stacked_rasters_n[[keep]]

#prefix
model_name <- "model4"

#preds_model4 <- mask(preds_model4, hull_vect)

# data
myBiomodData <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(coords)),
  resp.xy = coords,
  resp.name = species_name,
  expl.var = preds_model4,
  PA.nb.rep = 1,
  PA.nb.absences = 10000,
  PA.strategy = 'random'
)

models_out[[model_name]] <- models_function(myBiomodData,model_name, preds_model4)

#################################################################
#### model 5
#################################################################
#
# Model 5 (bio3 + clc26 + clc”else” – bio4)
#1.	Slope
#2.	Dem
#3.	Bio3
#4.	Bio8
#5.	Bio9
#6.	Bio13
#7.	Bio14
#8.	Percentage of clc 26
#9.	Sum percentage of clc (23,24,25,27,29,31,32)

keep <- c(
          "eudem_dem_3035_europe",
          "eudem_slop_3035_europe",
          "wc2.1_30s_bio_3",
          "wc2.1_30s_bio_8",
          "wc2.1_30s_bio_9",
          "wc2.1_30s_bio_13",
          "wc2.1_30s_bio_14",
          "PercFrac_26",
          "PercSuitClassElse"
          #"PercSuitClassAll"
)

preds_model5 <- num_stack[[keep]]
preds <- stacked_rasters_n[[keep]]

#prefix
model_name <- "model5"

#preds_model5 <- mask(preds_model5, hull_vect)

#data
myBiomodData <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(coords)),
  resp.xy = coords,
  resp.name = species_name,
  expl.var = preds_model5,
  PA.nb.rep = 1,
  PA.nb.absences = 10000,
  PA.strategy = 'random'
)

models_out[[model_name]] <- models_function(myBiomodData,model_name, preds=preds_model5)

#################################################################
#### model 6
#################################################################
#
#Model 6 (bio2 + clc26 + clc”else” – bio4)
#1.	Slope
#2.	Dem
#3.	Bio2
#4.	Bio8
#5.	Bio9
#6.	Bio13
#7.	Bio14
#8.	Percentage of clc 26
#9.	Sum percentage of clc (23,24,25,27,29,31,32)

keep <- c(
          "eudem_dem_3035_europe",
          "eudem_slop_3035_europe",
          "wc2.1_30s_bio_2",
          "wc2.1_30s_bio_8",
          "wc2.1_30s_bio_9",
          "wc2.1_30s_bio_13",
          "wc2.1_30s_bio_14",
          "PercFrac_26",
          "PercSuitClassElse"
          #"PercSuitClassAll"
)

preds_model6 <- num_stack[[keep]]
preds <- stacked_rasters_n[[keep]]

#prefix
#model_name <- "model6_RF"
model_name <- "model6_maxnet"

#preds_model6 <- mask(preds_model6, hull_vect)

# data
myBiomodData <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(coords)),
  resp.xy = coords,
  resp.name = species_name,
  expl.var = preds_model6,
  PA.nb.rep = 1,
  PA.nb.absences = 10000,
  PA.strategy = 'random'
)

models_out[[model_name]] <- models_function(myBiomodData,model_name, preds_model6)

#
#################################################################
#### model 7
#################################################################
#
#Model 7 (bio2 + clc26 + clc”else” – bio4)
#1.	Slope
#2.	Dem
#3.	Bio2
#4.	Bio8
#5.	Bio9
#6.	Bio13
#7.	Bio14
#8.	Percentage of clc 26
#9.	Sum percentage of clc (23,24,25,27,29,31,32)

keep <- c(
          "eudem_dem_3035_europe",
          "eudem_slop_3035_europe",
          "wc2.1_30s_bio_2",
          "wc2.1_30s_bio_8",
          "wc2.1_30s_bio_9",
          "wc2.1_30s_bio_13",
          "wc2.1_30s_bio_14",
          "PercFrac_26",
          "PercSuitClassElse"
          #"PercSuitClassAll"
)

preds_model7 <- num_stack[[keep]]
preds <- stacked_rasters_n[[keep]]

#prefix
model_name <- "model7"


parnassius_dist_b <- parnassius_dist |>
    st_make_valid() |>
    st_buffer(30000) |>
    st_union()


p_d <- vect(parnassius_dist_b)

preds_model7_c <- crop(preds_model7, p_d)

preds_model7_m <- mask(preds_model7_c, p_d)


# -------------------------
## check for NA before running the model
# -------------------------
which(is.na(terra::extract(preds_model7_m, p_apollo_points)))


png("../figures/rasters_model7.png",
    width = 6000,
    height = 1500,
    res=300,
    units="px")
par(mfrow = c(1,4),        # 1 row, 4 columns
    oma = c(1,1,1,1),      # outer margins
    mar = c(1,1,1,1),      # inner margins
    cex = 1.2)             # text size
plot(vect(p_apollo_points), legend = FALSE, main = "Points")
plot(preds_model7_m[[1]], legend = FALSE, main = "Raster")
plot(vect(parnassius_dist), legend = FALSE, main = "Action plan")
plot(vect(parnassius_dist_b), legend = FALSE, main = "Action plan buffer")
dev.off()


# data
myBiomodData <- BIOMOD_FormatingData(
  resp.var = rep(1, nrow(coords)),
  resp.xy = coords,
  resp.name = species_name,
  expl.var = preds_model7,
  PA.nb.rep = 1,
  PA.nb.absences = 10000,
  PA.dist.min = 10000,
  PA.strategy = 'disk'
)

models_out[[model_name]] <- models_function(myBiomodData,model_name, preds_model7_m)
#

# --------------------------------------------------
# save models to a list
# --------------------------------------------------
save(models_out, file = "../results/p_apollo_models_out.RData")

# --------------------------------------------------
# Favourite Reference Values
# --------------------------------------------------

# --------------------------------------------------
# MAXENT
# --------------------------------------------------

#library(MaxentVariableSelection)
#
#backgroundlocations <- system.file("extdata",
#"Backgrounddata.csv",
#package="MaxentVariableSelection")
#backgroundlocations <- read.csv(backgroundlocations,header=TRUE)
#head(backgroundlocations)
#
#
#VariableSelection( maxent="maxent.jar",
#outdir="OutputDirectory",
#gridfolder="BioORACLEVariables",
#occurrencelocations=system.file("extdata", "Occurrencedata.csv", package="MaxentVariableSelection"),
#backgroundlocations=system.file("extdata", "Backgrounddata.csv", package="MaxentVariableSelection"),
#additionalargs="nolinear noquadratic noproduct nothreshold noautofeature",
#contributionthreshold=5,
#correlationthreshold=0.7,
#betamultiplier=seq(2,6,0.5)
#)

# ------------------------------------
# Code graveyard
# ------------------------------------

