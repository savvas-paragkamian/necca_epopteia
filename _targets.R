library(targets)
library(tarchetypes)
library(yaml)

tar_option_set(
  packages = c(
    "dplyr",
    "readr",
    "readxl",
    "sf",
    "terra",
    "tibble",
    "tidyr",
    "stringr",
    "purrr",
    "units",
    "ggplot2",
    "ggnewscale"
  )
)

source("R/extract_occurrences.R")
source("R/extract_spatial.R")
source("R/transform.R")
source("R/load_maps.R")
source("R/load_official_outputs.R")
source("R/qc.R")
source("R/helper_functions.R")

a17 <- yaml::read_yaml("config/params.yml")

list(
  tar_target(a17_config, a17),
  
  tar_target(

    spatial_layers,
    read_spatial_layers(a17_config$inputs)
  ),
  
  tar_target(

    species_taxonomy,

    read_species_taxonomy(a17_config$inputs$taxonomy_curated)

  ),

  tar_target(

    gbif_occurrences,

    read_gbif_occurrences(

      path = a17_config$inputs$gbif_occurrences,

      greece_regions = spatial_layers$greece_regions

    )

  ),

  tar_target(

    e1x_mdpp_occurrences,

    read_e1x_mdpp_occurrences(a17_config$inputs$e1x_mdpp)

  ),

  tar_target(

    e1x_db_reference_occurrences,

    read_e1x_db_reference_occurrences(a17_config$inputs$e1x_ref)

  ),

  tar_target(

    e1x_db_sampling_occurrences,

    read_e1x_db_sampling_occurrences(a17_config$inputs$e1x_db)

  ),

  tar_target(

    e2x_occurrences,

    read_e2x_occurrences(a17_config$inputs$e2x_db)

  ),

  tar_target(

    private_occurrences,

    read_private_occurrences(a17_config$inputs$private_occurrences)

  ),

  tar_target(

    e2x_ref_unio_crassus,

    read_unio_crassus_occurrences(a17_config$inputs$e2x_ref_unio_crassus)
 

  ),

  tar_target(

    e2x_ref_stenobothrus_eurasius,

    read_stenobothrus_occurrences(a17_config$inputs$e2x_ref_stenobothrus_eurasius)

  ),

  tar_target(

    necca_redlist_points_occurrences,

    read_necca_redlist_points_occurrences(a17_config$inputs$necca_redlist_points)

  ),
  
  tar_target(
    national_report_distribution_occurrences,
    read_national_report_distribution_occurrences(
      distribution_path =
        a17_config$inputs$national_report_distribution,
  
      sensitive_distribution_path =
        a17_config$inputs$national_report_distribution_sensitive,
  
      species_taxonomy = species_taxonomy
    )
  ),

  tar_target(
    p_apollo_action_plan_occurrences,
    read_p_apollo_action_plan_occurrences(
      path = a17_config$inputs$p_apollo_action_plan,
      eea_grid_10km = spatial_layers$eea_grid_10km
    )
  ),

  # --- Transform ---

  tar_target(
    national_report_distribution_grid,
    build_national_report_distribution_grid_occurrences(
      species_distribution = national_report_distribution_occurrences,
      eea_grid_10km        = spatial_layers$eea_grid_10km
    )
  ),

  tar_target(
    species_occurrences_invertebrates,
    combine_all_occurrences(
      gbif_occurrences              = gbif_occurrences,
      e1x_mdpp_occurrences          = e1x_mdpp_occurrences,
      e1x_db_sampling_occurrences   = e1x_db_sampling_occurrences,
      e1x_db_reference_occurrences  = e1x_db_reference_occurrences,
      e2x_occurrences               = e2x_occurrences,
      necca_redlist_points_occurrences = necca_redlist_points_occurrences,
      e2x_ref_unio_crassus          = e2x_ref_unio_crassus,
      e2x_ref_stenobothrus_eurasius = e2x_ref_stenobothrus_eurasius,
      private_occurrences           = private_occurrences
    )
  ),

  tar_target(
    species_samples_art17,
    filter_art17_occurrences(
      species_occurrences_invertebrates = species_occurrences_invertebrates,
      species_taxonomy                  = species_taxonomy
    )
  ),

  tar_target(
    species_samples_eea,
    assign_eea_grid_10km(
      species_samples_art17 = species_samples_art17,
      eea_grid_10km         = spatial_layers$eea_grid_10km
    )
  ),

  tar_target(
    species_samples_presence_minimum,
    build_presence_minimum(
      species_samples_eea              = species_samples_eea,
      p_apollo_action_plan_occurrences = p_apollo_action_plan_occurrences,
      national_report_distribution_grid = national_report_distribution_grid
    )
  ),

  tar_target(
    species_samples_presence_elevation,
    enrich_with_elevation(
      species_samples_presence_minimum = species_samples_presence_minimum,
      eu_dem_path                      = a17_config$inputs$eu_dem
    )
  ),

  tar_target(
    species_samples_presence_border,
    enrich_with_border_distance(
      species_samples_presence_sf = species_samples_presence_elevation,
      greece_regions              = spatial_layers$greece_regions
    )
  ),

  tar_target(
    species_samples_presence_dist_flags,
    apply_distribution_filters(
      species_samples_presence_etrs89 = species_samples_presence_border
    )
  ),

  tar_target(
    species_samples_presence_pop,
    apply_population_filters(
      species_samples_presence_dist = species_samples_presence_dist_flags,
      eea_grid_1km                  = spatial_layers$eea_grid_1km
    )
  ),

  tar_target(
    species_samples_presence_final,
    build_presence_final(
      species_samples_presence_pop = species_samples_presence_pop
    )
  ),

  tar_target(
    species_samples_presence_final_private,
    build_presence_final_private(
      species_samples_presence_pop = species_samples_presence_pop
    )
  ),

  # --- Load (Maps) ---

  tar_target(
    map_borders_quality,
    save_borders_quality_map(
      species_samples_presence_border = species_samples_presence_border,
      greece_regions                  = spatial_layers$greece_regions,
      eea_grid_10km                   = spatial_layers$eea_grid_10km,
      path = file.path(a17_config$paths$maps_dir,
                       "map_occurrences_borders_filtering.png")
    ),
    format = "file"
  ),

  tar_target(
    map_natura_base,
    save_natura_base_map(
      greece_regions = spatial_layers$greece_regions,
      natura2000     = spatial_layers$natura2000,
      path = file.path(a17_config$paths$maps_dir, "map_natura.png")
    ),
    format = "file"
  ),

  tar_target(
    map_art17_overview,
    save_art17_overview_map(
      species_samples_presence_pop = species_samples_presence_pop,
      greece_regions               = spatial_layers$greece_regions,
      natura2000                   = spatial_layers$natura2000,
      path = file.path(a17_config$paths$maps_dir,
                       "map_art17_invertebrates_natura.png")
    ),
    format = "file"
  ),

  tar_target(
    maps_species_occurrences,
    save_species_occurrence_maps(
      species_samples_presence_pop = species_samples_presence_pop,
      greece_regions               = spatial_layers$greece_regions,
      natura2000                   = spatial_layers$natura2000,
      maps_dir                     = a17_config$paths$maps_dir
    ),
    format = "file"
  ),

  tar_target(
    species_range,
    compute_species_range(
      species_samples_presence_pop = species_samples_presence_pop,
      eea_grid_10km                = spatial_layers$eea_grid_10km
    )
  ),

  tar_target(
    maps_species_range,
    save_species_range_maps(
      species_range                = species_range,
      species_samples_presence_pop = species_samples_presence_pop,
      greece_regions               = spatial_layers$greece_regions,
      natura2000                   = spatial_layers$natura2000,
      maps_dir                     = a17_config$paths$maps_dir
    ),
    format = "file"
  ),

  tar_target(
    file_species_range_shp,
    save_species_range_shp(
      species_range = species_range,
      path = file.path(a17_config$paths$results_dir,
                       "species_range", "species_range.shp")
    ),
    format = "file"
  ),

  tar_target(
    maps_species_distribution,
    save_species_distribution_maps(
      species_samples_presence_pop = species_samples_presence_pop,
      eea_grid_10km                = spatial_layers$eea_grid_10km,
      greece_regions               = spatial_layers$greece_regions,
      natura2000                   = spatial_layers$natura2000,
      maps_dir                     = a17_config$paths$maps_dir
    ),
    format = "file"
  ),

  # --- Load (Official Outputs) ---

  tar_target(
    file_occurrences_invertebrates,
    save_occurrences_tsv(
      species_occurrences_invertebrates = species_occurrences_invertebrates,
      path = a17_config$outputs$species_occurrences_invertebrates
    ),
    format = "file"
  ),

  tar_target(
    file_art17_all,
    save_art17_tsv(
      species_samples_art17 = species_samples_art17,
      path = a17_config$outputs$species_samples_art17_all
    ),
    format = "file"
  ),

  tar_target(
    file_presence_final,
    save_presence_final_tsv(
      species_samples_presence_final = species_samples_presence_final,
      path = a17_config$outputs$species_samples_presence_final
    ),
    format = "file"
  ),

  tar_target(
    file_distributions_final,
    save_distributions_tsv(
      species_samples_presence_final = species_samples_presence_final,
      species_range                  = species_range,
      path = a17_config$outputs$distributions_presence_final
    ),
    format = "file"
  ),

  tar_target(
    file_populations_final,
    save_populations_tsv(
      species_samples_presence_final = species_samples_presence_final,
      path = a17_config$outputs$populations_presence_final
    ),
    format = "file"
  )

)



