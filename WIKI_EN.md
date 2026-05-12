# User Guide — necca_epopteia_pipeline

Detailed manual for using and extending the invertebrate data processing pipeline
for the Article 17 reporting (Directive 92/43/EEC), period 2019–2024.

---

## Contents

1. [Prerequisites](#1-prerequisites)
2. [Installation](#2-installation)
3. [One-time data preparation](#3-one-time-data-preparation)
4. [Folder structure](#4-folder-structure)
5. [Configuration — `config/params.yml`](#5-configuration--configparamsyml)
6. [Running the pipeline](#6-running-the-pipeline)
7. [Inspecting results](#7-inspecting-results)
8. [How to add a new data source](#8-how-to-add-a-new-data-source)
9. [Common errors](#9-common-errors)

---

## 1. Prerequisites

Before you begin, make sure the following are installed:

| Software | Version | Purpose |
|----------|---------|---------|
| R | ≥ 4.5.3 | Programming language |
| `renv` (R package) | any | Package management |
| Podman or Docker | any | Container (optional, for full reproducibility) |
| Git | any | Cloning the repository |

---

## 2. Installation

### 2.1 Clone the repository

```bash
git clone https://github.com/savvas-paragkamian/necca_epopteia.git
cd necca_epopteia
```

### 2.2 Restore the R environment

Open R inside the project folder and run:

```r
renv::restore()
```

This command automatically installs all packages at the exact versions defined
in `renv.lock`. It requires an internet connection and may take several minutes
the first time.

> **Note:** The `.Rprofile` file at the project root automatically activates
> `renv` every time you open R in this folder.

### 2.3 Alternative — Container (recommended)

For full reproducibility independent of the operating system:

```bash
podman build -t myproj-rgeo:4.5.3 -f .devcontainer/Containerfile .

podman run --rm -it \
    --userns=keep-id \
    -v "$PWD":/workspaces/project:Z \
    myproj-rgeo:4.5.3
```

Inside the container, restore the packages:

```r
renv::restore()
```

If the machine does not have enough memory:

```bash
podman machine stop
podman machine set --memory 12288 --cpus 5
podman machine start
```

### 2.4 Placing data files

Data files are **not** included in the repository. They must be placed manually
in the folders defined in `config/params.yml`:

- Occurrence data → `data/raw/`
- Spatial layers → `data/spatial/`

---

## 3. One-time data preparation

Before running the pipeline for the first time, three preparation steps must be
completed interactively. All helper functions live in `R/helper_functions.R`.

### 3.1 Taxonomic name verification

Species names must be verified against authoritative databases and **manually
curated** to produce `data/raw/species_taxonomy_curated.tsv`, which is the
taxonomy input for the pipeline. This step cannot be automated end-to-end —
the biologist must review the matches and decide on accepted canonical names.

```r
source("R/helper_functions.R")

verify_species_taxonomy(
  species_names = species_names_combined,
  output_path   = "results/gnr_species_verifier.tsv"
)
```

Names are checked against Catalogue of Life (1), WoRMS (9), GBIF (11) and
EOL (12) with `all_matches = TRUE`. After reviewing the output TSV, update
`data/raw/species_taxonomy_curated.tsv` with the accepted names.

### 3.2 GBIF occurrence download

GBIF credentials must be stored in `~/.Renviron`:

```
GBIF_USER=your_username
GBIF_PWD=your_password
GBIF_EMAIL=your_email
```

Then run interactively:

```r
source("R/helper_functions.R")

keys <- get_gbif_taxon_keys(species_names_combined)
key  <- request_gbif_download(keys, country = "GR")
import_gbif_download(key, output_path = "data/raw/gbif_invertebrate_species_occ.tsv")
```

Or as a single call:

```r
download_gbif_occurrences(
  species_names = species_names_combined,
  output_path   = "data/raw/gbif_invertebrate_species_occ.tsv"
)
```

### 3.3 Cropping large rasters to Greece

The full EU DEM mosaic (~20 GB, EPSG:3035) must be cropped to Greece before
the pipeline can use it:

```r
source("R/helper_functions.R")

greece_regions <- sf::st_read("data/spatial/gadm41_GRC_shp/gadm41_GRC_2.shp")

crop_eu_dem_to_greece(
  eu_dem_path    = "/path/to/full/eudem_dem_3035_europe.tif",
  greece_regions = greece_regions,
  output_path    = "data/spatial/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_3035_europe.tif"
)
```

For any other large raster (WorldClim, CORINE, etc.):

```r
crop_raster_to_extent(
  raster_path = "/path/to/source.tif",
  extent_sf   = greece_regions,
  output_path = "data/spatial/output/cropped.tif",
  target_crs  = 3035   # NULL to keep source CRS
)
```

---

## 4. Folder structure


```
necca_epopteia_pipeline/
├── _targets.R                  # Pipeline definition (entry point)
├── config/
│   └── params.yml              # All file paths and parameters
├── R/
│   ├── extract_occurrences.R   # Reading occurrence data
│   ├── extract_spatial.R       # Reading spatial layers
│   ├── transform.R             # Filtering, enrichment, flags
│   ├── load_maps.R             # Map generation (PNG)
│   ├── load_official_outputs.R # TSV report generation
│   ├── helper_functions.R      # Spatial utility functions
│   └── qc.R                    # Quality control (under development)
├── data/
│   ├── raw/                    # Input files (not tracked by git)
│   ├── spatial/                # Spatial reference layers
│   └── derived/                # Intermediate results
├── results/
│   ├── official/               # Final TSVs for the Article 17 report
│   ├── maps/                   # PNG maps
│   └── species_range/          # Range distribution shapefile
├── renv.lock                   # Pinned package versions
└── WIKI_EN.md                  # This file
```

---

## 5. Configuration — `config/params.yml`

The file `config/params.yml` is the **only place** where file paths are declared.
There are no hardcoded paths anywhere in the R source code.

File structure:

```yaml
reporting:
  framework: "A17"
  period: "2019-2024"

paths:
  raw_data_dir: "data/raw"
  results_dir: "results"
  maps_dir: "results/maps"

inputs:
  taxonomy_curated: "data/raw/species_taxonomy_curated.tsv"
  gbif_occurrences: "data/raw/gbif_invertebrate_species_occ.tsv"
  e1x_mdpp: "data/raw/Ε1Χ_ΒΔ_ΠΡΩΤΟΓΕΝΩΝ_ΦΔ+ΜΔΠΠ_2014-2024_v8.xlsx"
  # ... remaining sources ...
  greece_regions: "data/spatial/gadm41_GRC_shp/gadm41_GRC_2.shp"
  eea_grid_10km:  "data/spatial/eea_reference_grid/gr_10km.shp"
  eu_dem: "data/spatial/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_3035_europe.tif"

outputs:
  species_occurrences_invertebrates: "data/derived/species_occurrences_invertebrates.tsv"
  species_samples_presence_final: "results/official/species_samples_presence_final.tsv"
  # ...
```

If your files are in a different location, change only the paths here — no other
code needs to be modified.

---

## 6. Running the pipeline

### 6.1 Full run

Open R in the project folder and run:

```r
targets::tar_make()
```

The pipeline automatically executes all Extract → Transform → Load steps in the
correct order. Each step is cached in the `_targets/` folder. If you re-run the
command without any changes to the input data, `targets` skips the already
computed steps.

### 6.2 Parallel execution

For faster execution on a multi-core system:

```r
library(future)
plan(multisession)
targets::tar_make(callr_function = NULL)
```

### 6.3 Partial execution

You can run only a specific step:

```r
targets::tar_make("species_samples_eea")
```

`targets` will automatically run any upstream steps that are needed.

### 6.4 Useful management commands

```r
# Which steps need to re-run?
targets::tar_outdated()

# Visualise the dependency graph in the browser
targets::tar_visnetwork()

# Read the result of a step without re-running it
targets::tar_read(species_samples_presence_final)

# Full reset — delete the cache
targets::tar_destroy()
```

---

## 7. Inspecting results

### 7.1 Final output files

After a successful run, output files are located at the paths defined under
`outputs:` in `config/params.yml`:

| File | Content |
|------|---------|
| `data/derived/species_occurrences_invertebrates.tsv` | All records from all sources combined |
| `data/derived/species_samples_art17_all.tsv` | Records filtered to Annex II species |
| `results/official/species_samples_presence_final.tsv` | Final presence dataset (no private records) |
| `results/official/distributions_presence_final.tsv` | Distribution per species and 10 km cell |
| `results/official/populations_presence_final.tsv` | Population per species and 1 km cell |
| `results/species_range/species_range.shp` | Species range polygons |
| `results/maps/` | Per-species and overview PNG maps |

### 7.2 Inspecting intermediate steps

Any intermediate result can be loaded directly into memory:

```r
# View records after EEA grid assignment
targets::tar_read(species_samples_eea)

# View records after distribution filters
targets::tar_read(species_samples_presence_dist_flags)

# Inspect the range computation
targets::tar_read(species_range)
```

---

## 8. How to add a new data source

Adding a new occurrence data source requires changes to **four files** in a
specific order. Follow the steps below.

---

### Step 1 — Add the file path to `config/params.yml`

Open `config/params.yml` and add a new line under `inputs:`:

```yaml
inputs:
  # ... existing sources ...
  my_new_source: "data/raw/my_new_data_file.xlsx"
```

The key (e.g. `my_new_source`) is how the file will be referenced in the code.
The path is relative to the project root.

---

### Step 2 — Write a reader function in `R/extract_occurrences.R`

Add a new function at the end of the file. The function must return a
`data.frame` with at least the Darwin Core columns:

```r
read_my_new_source_occurrences <- function(path) {
  readxl::read_xlsx(path) |>
    dplyr::rename(
      submittedName    = ScientificName,   # adapt to your column names
      decimalLatitude  = Latitude,
      decimalLongitude = Longitude
    ) |>
    dplyr::mutate(
      datasetName     = "My_New_Source",     # unique source identifier
      basisOfRecord   = "MATERIAL_SAMPLE",   # or MaterialCitation, HUMAN_OBSERVATION
      recordNumber    = as.character(ID),
      collectionCode  = basename(path),
      individualCount = as.numeric(Count)
    )
}
```

**Required output columns** (Darwin Core):

| Column | Type | Description |
|--------|------|-------------|
| `submittedName` | character | Scientific name as it appears in the source |
| `decimalLatitude` | numeric | Latitude in WGS84 |
| `decimalLongitude` | numeric | Longitude in WGS84 |
| `datasetName` | character | Unique source identifier |
| `recordNumber` | character | Unique record identifier |
| `collectionCode` | character | Source filename |
| `basisOfRecord` | character | Record type |
| `individualCount` | numeric | Number of individuals (NA if unknown) |

---

### Step 3 — Add a target to `_targets.R`

Open `_targets.R` and add a new `tar_target()` in the Extract section,
alongside the other `read_*_occurrences` targets:

```r
tar_target(
  my_new_source_occurrences,
  read_my_new_source_occurrences(a17_config$inputs$my_new_source)
),
```

The target name (`my_new_source_occurrences`) will be used in the next step.

---

### Step 4 — Integrate into `combine_all_occurrences()`

The new source must be wired in at **two points**:

#### 4a. In `_targets.R` — add the argument to the call

Find the `species_occurrences_invertebrates` target and add the new source:

```r
tar_target(
  species_occurrences_invertebrates,
  combine_all_occurrences(
    gbif_occurrences              = gbif_occurrences,
    e1x_mdpp_occurrences          = e1x_mdpp_occurrences,
    # ... remaining sources ...
    my_new_source_occurrences     = my_new_source_occurrences   # ← new line
  )
),
```

#### 4b. In `R/transform.R` — add the parameter to the function

Find the `combine_all_occurrences()` function and:

1. Add the new argument to the parameter list:

```r
combine_all_occurrences <- function(
  gbif_occurrences,
  e1x_mdpp_occurrences,
  # ... remaining parameters ...
  my_new_source_occurrences       # ← new parameter
) {
```

2. Add the new object to the `list()` inside the function:

```r
  list(
    gbif_occurrences,
    e1x_mdpp_occurrences,
    # ... remaining ...
    my_new_source_occurrences     # ← new line
  ) |>
    purrr::map(~ dplyr::select(.x, dplyr::all_of(cols))) |>
    dplyr::bind_rows()
```

---

### Step 5 — Verify (optional but recommended)

Run only the new targets for a quick check before a full pipeline run:

```r
# Test reading the new source
targets::tar_make("my_new_source_occurrences")
targets::tar_read("my_new_source_occurrences")

# Test combining all sources
targets::tar_make("species_occurrences_invertebrates")
targets::tar_read("species_occurrences_invertebrates") |> dplyr::count(datasetName)
```

Confirm that the new source appears in the `datasetName` column of the result.

---

### Step 6 — Custom filtering rules (if needed)

If the new source requires special treatment in the filtering stage (e.g.
exclusion from population counts, or different flag logic), open `R/transform.R`
and modify the appropriate functions:

- **`apply_distribution_filters()`** — for rules affecting the `includeDistribution` column
- **`apply_population_filters()`** — for rules affecting the `includePopulation` column

Example — exclude the new source from population counts:

```r
apply_population_filters <- function(species_samples_presence_dist, eea_grid_1km) {
  # ... existing logic ...
  dplyr::mutate(includePopulation = dplyr::if_else(
    datasetName == "My_New_Source",
    FALSE, includePopulation
  ))
}
```

---

### Summary of changes for a new data source

| File | Change |
|------|--------|
| `config/params.yml` | New line under `inputs:` with the file path |
| `R/extract_occurrences.R` | New `read_*()` function returning Darwin Core columns |
| `_targets.R` (Extract) | New `tar_target()` calling the reader function |
| `_targets.R` (Transform) | New argument added to the `combine_all_occurrences()` call |
| `R/transform.R` | New parameter and object added to `combine_all_occurrences()` |
| `R/transform.R` (if needed) | Filtering rules in `apply_*_filters()` |

---

## 9. Common errors

### `Error in loadNamespace: there is no package called 'targets'`

Packages have not been installed. Run:

```r
renv::restore()
```

### `external pointer is not valid`

This occurs if a `terra` raster object has been stored as a cached target.
In this pipeline, the EU DEM is **never** stored as a target object — only its
file path is passed. Do not call `terra::rast()` inside a `tar_target()` — call
it inside the function that needs the raster.

### Missing column from new source (`Error: Column X not found`)

The new `read_*()` function is not returning all required Darwin Core columns.
Check the table in Step 2 and make sure all eight columns are present before the
function returns.

### `tar_make()` does not re-run a step that was changed

`targets` tracks code and data via hashing. If you changed an R file but the
step does not appear as outdated, check with:

```r
targets::tar_outdated()
targets::tar_invalidate("my_target_name")  # force re-run
```
