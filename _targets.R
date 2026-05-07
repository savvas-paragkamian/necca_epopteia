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
    "units"
  )
)

source("R/extract_occurrences.R")
source("R/extract_spatial.R")
source("R/transform_enrichment.R")
source("R/transform_spatial.R")
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
        eea_grid_10km = reference_layers$eea_grid_10km
      )
    )

)



