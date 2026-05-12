# Οδηγός Χρήσης — necca_epopteia_pipeline

Αναλυτικό εγχειρίδιο για τη χρήση και επέκταση του αγωγού επεξεργασίας δεδομένων
ασπόνδυλων για την αναφορά Άρθρου 17 (Οδηγία 92/43/ΕΟΚ), περίοδος 2019–2024.

---

## Περιεχόμενα

1. [Προαπαιτούμενα](#1-προαπαιτούμενα)
2. [Εγκατάσταση](#2-εγκατάσταση)
3. [Προετοιμασία δεδομένων μίας-χρήσης](#3-προετοιμασία-δεδομένων-μίας-χρήσης)
4. [Δομή φακέλων](#4-δομή-φακέλων)
5. [Ρύθμιση — `config/params.yml`](#5-ρύθμιση--configparamsyml)
6. [Εκτέλεση του αγωγού](#6-εκτέλεση-του-αγωγού)
7. [Επιθεώρηση αποτελεσμάτων](#7-επιθεώρηση-αποτελεσμάτων)
8. [Πώς να προσθέσετε νέα πηγή δεδομένων](#8-πώς-να-προσθέσετε-νέα-πηγή-δεδομένων)
9. [Συνηθισμένα σφάλματα](#9-συνηθισμένα-σφάλματα)

---

## 1. Προαπαιτούμενα

Πριν ξεκινήσετε, βεβαιωθείτε ότι έχετε εγκατεστημένα:

| Λογισμικό | Έκδοση | Σκοπός |
|-----------|--------|--------|
| R | ≥ 4.5.3 | Γλώσσα προγραμματισμού |
| `renv` (R πακέτο) | οποιαδήποτε | Διαχείριση πακέτων |
| Podman ή Docker | οποιαδήποτε | Container (προαιρετικό, για πλήρη αναπαραγωγιμότητα) |
| Git | οποιαδήποτε | Κλωνοποίηση αποθετηρίου |

---

## 2. Εγκατάσταση

### 2.1 Κλωνοποίηση αποθετηρίου

```bash
git clone https://github.com/savvas-paragkamian/necca_epopteia.git
cd necca_epopteia
```

### 2.2 Αποκατάσταση περιβάλλοντος R

Ανοίξτε R μέσα στον φάκελο του project και εκτελέστε:

```r
renv::restore()
```

Αυτή η εντολή εγκαθιστά αυτόματα όλα τα πακέτα στις ακριβείς εκδόσεις που
ορίζονται στο `renv.lock`. Χρειάζεται σύνδεση στο διαδίκτυο και μπορεί να
διαρκέσει αρκετά λεπτά την πρώτη φορά.

> **Σημείωση:** Το αρχείο `.Rprofile` στη ρίζα του project ενεργοποιεί
> αυτόματα το `renv` κάθε φορά που ανοίγετε R σε αυτόν τον φάκελο.

### 2.3 Εναλλακτικά — Container (συνιστάται)

Για πλήρη αναπαραγωγιμότητα ανεξάρτητη από το λειτουργικό σύστημα:

```bash
podman build -t myproj-rgeo:4.5.3 -f .devcontainer/Containerfile .

podman run --rm -it \
    --userns=keep-id \
    -v "$PWD":/workspaces/project:Z \
    myproj-rgeo:4.5.3
```

Μέσα στο container, αποκαταστήστε τα πακέτα:

```r
renv::restore()
```

Αν ο υπολογιστής δεν έχει αρκετή μνήμη:

```bash
podman machine stop
podman machine set --memory 12288 --cpus 5
podman machine start
```

### 2.4 Τοποθέτηση δεδομένων

Τα αρχεία δεδομένων **δεν** συμπεριλαμβάνονται στο αποθετήριο.
Πρέπει να τοποθετηθούν χειροκίνητα στους φακέλους που ορίζονται στο
`config/params.yml`:

- Δεδομένα εμφανίσεων → `data/raw/`
- Χωρικά επίπεδα → `data/spatial/`

---

## 3. Προετοιμασία δεδομένων μίας-χρήσης

Πριν την πρώτη εκτέλεση του αγωγού απαιτούνται τρία προπαρασκευαστικά βήματα
που εκτελούνται διαδραστικά. Όλες οι βοηθητικές συναρτήσεις βρίσκονται στο
`R/helper_functions.R`.

### 3.1 Επαλήθευση ταξινόμησης ειδών

Τα ονόματα των ειδών πρέπει να επαληθευτούν έναντι αυθεντικών βάσεων
δεδομένων και να υποστούν **χειροκίνητη επιμέλεια** για να παραχθεί το
`data/raw/species_taxonomy_curated.tsv`, το οποίο αποτελεί είσοδο ταξινόμησης
για τον αγωγό. Το βήμα αυτό δεν μπορεί να αυτοματοποιηθεί πλήρως — ο
βιολόγος πρέπει να επιθεωρήσει τις αντιστοιχίσεις και να αποφασίσει για τα
αποδεκτά canonical ονόματα.

```r
source("R/helper_functions.R")

verify_species_taxonomy(
  species_names = species_names_combined,
  output_path   = "results/gnr_species_verifier.tsv"
)
```

Τα ονόματα ελέγχονται έναντι Catalogue of Life (1), WoRMS (9), GBIF (11) και
EOL (12) με `all_matches = TRUE`. Μετά την επιθεώρηση του αρχείου εξόδου,
ενημερώστε το `data/raw/species_taxonomy_curated.tsv` με τα αποδεκτά ονόματα.

### 3.2 Λήψη καταγραφών GBIF

Τα διαπιστευτήρια GBIF πρέπει να αποθηκευτούν στο `~/.Renviron`:

```
GBIF_USER=your_username
GBIF_PWD=your_password
GBIF_EMAIL=your_email
```

Στη συνέχεια εκτελέστε διαδραστικά:

```r
source("R/helper_functions.R")

keys <- get_gbif_taxon_keys(species_names_combined)
key  <- request_gbif_download(keys, country = "GR")
import_gbif_download(key, output_path = "data/raw/gbif_invertebrate_species_occ.tsv")
```

Ή ως ενιαία κλήση:

```r
download_gbif_occurrences(
  species_names = species_names_combined,
  output_path   = "data/raw/gbif_invertebrate_species_occ.tsv"
)
```

### 3.3 Περικοπή μεγάλων rasters στην Ελλάδα

Το πλήρες EU DEM mosaic (~20 GB, EPSG:3035) πρέπει να περικοπεί στα όρια
της Ελλάδας πριν τη χρήση του από τον αγωγό:

```r
source("R/helper_functions.R")

greece_regions <- sf::st_read("data/spatial/gadm41_GRC_shp/gadm41_GRC_2.shp")

crop_eu_dem_to_greece(
  eu_dem_path    = "/path/to/full/eudem_dem_3035_europe.tif",
  greece_regions = greece_regions,
  output_path    = "data/spatial/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_3035_europe.tif"
)
```

Για οποιοδήποτε άλλο μεγάλο raster (WorldClim, CORINE κτλ.):

```r
crop_raster_to_extent(
  raster_path = "/path/to/source.tif",
  extent_sf   = greece_regions,
  output_path = "data/spatial/output/cropped.tif",
  target_crs  = 3035   # NULL για διατήρηση του αρχικού CRS
)
```

---

## 4. Δομή φακέλων

```
necca_epopteia_pipeline/
├── _targets.R                  # Ορισμός αγωγού (σημείο εκκίνησης)
├── config/
│   └── params.yml              # Όλες οι διαδρομές αρχείων και παράμετροι
├── R/
│   ├── extract_occurrences.R   # Ανάγνωση δεδομένων παρουσίας
│   ├── extract_spatial.R       # Ανάγνωση χωρικών επιπέδων
│   ├── transform.R             # Φιλτράρισμα, εμπλουτισμός, σημαίες
│   ├── load_maps.R             # Παραγωγή χαρτών (PNG)
│   ├── load_official_outputs.R # Παραγωγή TSV αναφορών
│   ├── helper_functions.R      # Βοηθητικές χωρικές συναρτήσεις
│   └── qc.R                    # Ποιοτικός έλεγχος (υπό ανάπτυξη)
├── data/
│   ├── raw/                    # Αρχεία εισόδου (δεν ανεβαίνουν στο git)
│   ├── spatial/                # Χωρικά επίπεδα αναφοράς
│   └── derived/                # Ενδιάμεσα αποτελέσματα
├── results/
│   ├── official/               # Τελικά TSV για την αναφορά Άρθρου 17
│   ├── maps/                   # Χάρτες PNG
│   └── species_range/          # Shapefile εύρους κατανομής
├── renv.lock                   # Καρφιτσωμένες εκδόσεις πακέτων
└── WIKI.md                     # Αυτό το αρχείο
```

---

## 5. Ρύθμιση — `config/params.yml`

Το αρχείο `config/params.yml` είναι ο **μοναδικός τόπος** όπου δηλώνονται
διαδρομές αρχείων. Δεν υπάρχουν hardcoded διαδρομές μέσα στον κώδικα R.

Η δομή του αρχείου:

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
  # ... υπόλοιπες πηγές ...
  greece_regions: "data/spatial/gadm41_GRC_shp/gadm41_GRC_2.shp"
  eea_grid_10km:  "data/spatial/eea_reference_grid/gr_10km.shp"
  eu_dem: "data/spatial/EU_DEM_mosaic_5deg_gr/crop_eudem_dem_3035_europe.tif"

outputs:
  species_occurrences_invertebrates: "data/derived/species_occurrences_invertebrates.tsv"
  species_samples_presence_final: "results/official/species_samples_presence_final.tsv"
  # ...
```

Αν τα αρχεία σας βρίσκονται σε διαφορετική τοποθεσία, αλλάξτε μόνο
τις διαδρομές εδώ — ο υπόλοιπος κώδικας δεν χρειάζεται καμία τροποποίηση.

---

## 6. Εκτέλεση του αγωγού

### 5.1 Πλήρης εκτέλεση

Ανοίξτε R στον φάκελο του project και εκτελέστε:

```r
targets::tar_make()
```

Ο αγωγός εκτελεί αυτόματα όλα τα βήματα Extract → Transform → Load με τη
σωστή σειρά. Κάθε βήμα αποθηκεύεται στον φάκελο `_targets/`. Αν
επανεκτελέσετε την εντολή χωρίς να έχουν αλλάξει τα δεδομένα εισόδου,
το `targets` παραλείπει τα ήδη υπολογισμένα βήματα.

### 5.2 Παράλληλη εκτέλεση

Για ταχύτερη εκτέλεση σε πολυπύρηνο σύστημα:

```r
library(future)
plan(multisession)
targets::tar_make(callr_function = NULL)
```

### 5.3 Μερική εκτέλεση

Μπορείτε να εκτελέσετε μόνο ένα συγκεκριμένο βήμα:

```r
targets::tar_make("species_samples_eea")
```

Το `targets` θα εκτελέσει αυτόματα και τα προηγούμενα βήματα αν χρειάζεται.

### 5.4 Χρήσιμες εντολές διαχείρισης

```r
# Ποια βήματα χρειάζονται επανεκτέλεση;
targets::tar_outdated()

# Οπτικοποίηση γράφου εξαρτήσεων στο browser
targets::tar_visnetwork()

# Ανάγνωση αποτελέσματος ενός βήματος χωρίς επανεκτέλεση
targets::tar_read(species_samples_presence_final)

# Πλήρης επαναφορά — διαγραφή της κρυφής μνήμης
targets::tar_destroy()
```

---

## 7. Επιθεώρηση αποτελεσμάτων

### 6.1 Τελικά αρχεία εξόδου

Μετά την επιτυχή εκτέλεση, τα αρχεία εξόδου βρίσκονται στις διαδρομές
που ορίζονται στο `config/params.yml` υπό `outputs:`:

| Αρχείο | Περιεχόμενο |
|--------|-------------|
| `data/derived/species_occurrences_invertebrates.tsv` | Σύνολο καταγραφών από όλες τις πηγές |
| `data/derived/species_samples_art17_all.tsv` | Καταγραφές φιλτραρισμένες στα είδη Παραρτήματος ΙΙ |
| `results/official/species_samples_presence_final.tsv` | Τελικό σύνολο παρουσιών (χωρίς ιδιωτικά) |
| `results/official/distributions_presence_final.tsv` | Κατανομή ανά είδος και κελί 10 km |
| `results/official/populations_presence_final.tsv` | Πληθυσμός ανά είδος και κελί 1 km |
| `results/species_range/species_range.shp` | Polygons εύρους κατανομής ειδών |
| `results/maps/` | Χάρτες PNG ανά είδος και συνολικοί |

### 6.2 Επιθεώρηση ενδιάμεσων βημάτων

Κάθε ενδιάμεσο αποτέλεσμα μπορεί να φορτωθεί απευθείας στη μνήμη:

```r
# Δείτε τις καταγραφές μετά το φιλτράρισμα EEA
targets::tar_read(species_samples_eea)

# Δείτε τις καταγραφές μετά τα φίλτρα κατανομής
targets::tar_read(species_samples_presence_dist_flags)

# Ελέγξτε τον υπολογισμό εύρους
targets::tar_read(species_range)
```

---

## 8. Πώς να προσθέσετε νέα πηγή δεδομένων

Η προσθήκη νέας πηγής δεδομένων εμφανίσεων απαιτεί αλλαγές σε **τέσσερα
αρχεία** με συγκεκριμένη σειρά. Ακολουθήστε τα παρακάτω βήματα.

---

### Βήμα 1 — Προσθήκη διαδρομής στο `config/params.yml`

Ανοίξτε το `config/params.yml` και προσθέστε νέα γραμμή στο τμήμα `inputs:`:

```yaml
inputs:
  # ... υπάρχουσες πηγές ...
  my_new_source: "data/raw/my_new_data_file.xlsx"
```

Το κλειδί (π.χ. `my_new_source`) θα χρησιμοποιηθεί για να αναφερθείτε στο
αρχείο μέσα από τον κώδικα. Η διαδρομή είναι σχετική ως προς τη ρίζα του
project.

---

### Βήμα 2 — Γραφή συνάρτησης ανάγνωσης στο `R/extract_occurrences.R`

Προσθέστε νέα συνάρτηση στο τέλος του αρχείου. Η συνάρτηση πρέπει να
επιστρέφει `data.frame` με τουλάχιστον τις στήλες Darwin Core:

```r
read_my_new_source_occurrences <- function(path) {
  # Διαβάστε το αρχείο σας
  readxl::read_xlsx(path) |>
    dplyr::rename(
      submittedName    = ScientificName,   # προσαρμόστε στις δικές σας στήλες
      decimalLatitude  = Latitude,
      decimalLongitude = Longitude
    ) |>
    dplyr::mutate(
      collectionCode    = "My_New_Source",      # μοναδικό όνομα πηγής
      basisOfRecord  = "MATERIAL_SAMPLE",    # ή MaterialCitation, HUMAN_OBSERVATION
      recordNumber   = as.character(ID),
      datasetName = basename(path),
      individualCount = as.numeric(Count)
    )
}
```

**Απαιτούμενες στήλες εξόδου** (Darwin Core):

| Στήλη | Τύπος | Περιγραφή |
|-------|-------|-----------|
| `submittedName` | character | Επιστημονικό όνομα όπως υπάρχει στην πηγή |
| `decimalLatitude` | numeric | Γεωγραφικό πλάτος WGS84 |
| `decimalLongitude` | numeric | Γεωγραφικό μήκος WGS84 |
| `collectionCode` | character | Μοναδικό όνομα πηγής |
| `recordNumber` | character | Μοναδικός αριθμός εγγραφής |
| `datasetName` | character | Όνομα αρχείου προέλευσης |
| `basisOfRecord` | character | Τύπος εγγραφής |
| `individualCount` | numeric | Αριθμός ατόμων (NA αν άγνωστο) |

---

### Βήμα 3 — Προσθήκη target στο `_targets.R`

Ανοίξτε το `_targets.R` και προσθέστε νέο `tar_target()` στο τμήμα
Extract, μαζί με τα υπόλοιπα `read_*_occurrences`:

```r
tar_target(
  my_new_source_occurrences,
  read_my_new_source_occurrences(a17_config$inputs$my_new_source)
),
```

Το όνομα του target (`my_new_source_occurrences`) θα χρησιμοποιηθεί στο
επόμενο βήμα.

---

### Βήμα 4 — Ενσωμάτωση στη `combine_all_occurrences()`

Η ενσωμάτωση της νέας πηγής γίνεται σε **δύο σημεία**:

#### 4α. Στο `_targets.R` — προσθήκη ορίσματος στην κλήση

Βρείτε το target `species_occurrences_invertebrates` και προσθέστε τη νέα
πηγή:

```r
tar_target(
  species_occurrences_invertebrates,
  combine_all_occurrences(
    gbif_occurrences              = gbif_occurrences,
    e1x_mdpp_occurrences          = e1x_mdpp_occurrences,
    # ... υπόλοιπες πηγές ...
    my_new_source_occurrences     = my_new_source_occurrences   # ← νέα γραμμή
  )
),
```

#### 4β. Στο `R/transform.R` — προσθήκη παραμέτρου στη συνάρτηση

Βρείτε τη συνάρτηση `combine_all_occurrences()` και:

1. Προσθέστε το νέο όρισμα στη λίστα παραμέτρων:

```r
combine_all_occurrences <- function(
  gbif_occurrences,
  e1x_mdpp_occurrences,
  # ... υπόλοιπες παράμετροι ...
  my_new_source_occurrences       # ← νέα παράμετρος
) {
```

2. Προσθέστε το νέο αντικείμενο στη λίστα `list()` εντός της συνάρτησης:

```r
  list(
    gbif_occurrences,
    e1x_mdpp_occurrences,
    # ... υπόλοιπα ...
    my_new_source_occurrences     # ← νέα γραμμή
  ) |>
    purrr::map(~ dplyr::select(.x, dplyr::all_of(cols))) |>
    dplyr::bind_rows()
```

---

### Βήμα 5 — Επαλήθευση (προαιρετικό αλλά συνιστάται)

Εκτελέστε μόνο τα νέα targets για γρήγορο έλεγχο πριν την πλήρη εκτέλεση:

```r
# Δοκιμή ανάγνωσης νέας πηγής
targets::tar_make("my_new_source_occurrences")
targets::tar_read("my_new_source_occurrences")

# Δοκιμή ενοποίησης
targets::tar_make("species_occurrences_invertebrates")
targets::tar_read("species_occurrences_invertebrates") |> dplyr::count(collectionCode)
```

Βεβαιωθείτε ότι η νέα πηγή εμφανίζεται στη στήλη `collectionCode` του
αποτελέσματος.

---

### Βήμα 6 — Ειδικοί κανόνες φιλτραρίσματος (αν χρειάζεται)

Αν η νέα πηγή χρειάζεται ειδική μεταχείριση στο φιλτράρισμα (π.χ.
αποκλεισμός από τον υπολογισμό πληθυσμού, ή διαφορετική λογική για τα
flags), ανοίξτε το `R/transform.R` και τροποποιήστε τις κατάλληλες
συναρτήσεις:

- **`apply_distribution_filters()`** — για κανόνες που αφορούν τη στήλη
  `includeDistribution`
- **`apply_population_filters()`** — για κανόνες που αφορούν τη στήλη
  `includePopulation`

Παράδειγμα — αποκλεισμός νέας πηγής από τον πληθυσμό:

```r
apply_population_filters <- function(species_samples_presence_dist, eea_grid_1km) {
  # ... υπάρχουσα λογική ...
  dplyr::mutate(includePopulation = dplyr::if_else(
    collectionCode == "My_New_Source",
    FALSE, includePopulation
  ))
}
```

---

### Σύνοψη αλλαγών για νέα πηγή δεδομένων

| Αρχείο | Αλλαγή |
|--------|--------|
| `config/params.yml` | Νέα γραμμή στο `inputs:` με τη διαδρομή του αρχείου |
| `R/extract_occurrences.R` | Νέα συνάρτηση `read_*()` που επιστρέφει Darwin Core στήλες |
| `_targets.R` (Extract) | Νέο `tar_target()` που καλεί τη συνάρτηση ανάγνωσης |
| `_targets.R` (Transform) | Προσθήκη νέου ορίσματος στην κλήση `combine_all_occurrences()` |
| `R/transform.R` | Προσθήκη νέας παραμέτρου και αντικειμένου στη `combine_all_occurrences()` |
| `R/transform.R` (αν χρειάζεται) | Κανόνες φιλτραρίσματος στις `apply_*_filters()` |

---

## 9. Συνηθισμένα σφάλματα

### `Error in loadNamespace: there is no package called 'targets'`

Τα πακέτα δεν έχουν εγκατασταθεί. Εκτελέστε:

```r
renv::restore()
```

### `external pointer is not valid`

Εμφανίζεται αν ένα αντικείμενο `terra` (raster) έχει αποθηκευτεί ως
cached target. Στον αγωγό αυτό, το EU DEM **δεν** αποθηκεύεται — περνιέται
μόνο η διαδρομή του αρχείου. Μην προσθέτετε `terra::rast()` μέσα σε
`tar_target()` — καλέστε το εντός της συνάρτησης που το χρειάζεται.

### Στήλη λείπει από νέα πηγή (`Error: Column X not found`)

Η νέα συνάρτηση `read_*()` δεν επιστρέφει όλες τις απαιτούμενες στήλες
Darwin Core. Ελέγξτε τη λίστα στο Βήμα 2 και βεβαιωθείτε ότι όλες οι
οκτώ στήλες υπάρχουν πριν επιστρέψει η συνάρτηση.

### `tar_make()` δεν τρέχει ξανά ένα βήμα που άλλαξε

Το `targets` βασίζεται στο hash του κώδικα και των δεδομένων. Αν αλλάξατε
το αρχείο R αλλά το βήμα δεν εμφανίζεται ως outdated, ελέγξτε με:

```r
targets::tar_outdated()
targets::tar_invalidate("my_target_name")  # αναγκαστική επανεκτέλεση
```
