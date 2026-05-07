# Αγωγός επεξεργασίας δεδομένων ασπόνδυλων — Ελλάδα, Άρθρο 17

Ροή εργασιών ανάλυσης και αναφοράς για την παρακολούθηση ασπόνδυλων στο πλαίσιο
της Οδηγίας για τους Οικοτόπους Natura 2000 (Άρθρο 17, Οδηγία 92/43/ΕΟΚ),
περίοδος αναφοράς 2019–2024. Χρηματοδοτείται από τον Οργανισμό Φυσικού
Περιβάλλοντος και Κλιματικής Αλλαγής (ΟΦΥΠΕΚΑ), Οκτώβριος 2024 – Ιούλιος 2025.

Ο αγωγός ενσωματώνει δεδομένα παρουσίας από 10 πηγές, εμπλουτίζει τις
καταγραφές με χωρικά επίπεδα (πλέγματα EEA, Natura 2000, Ψηφιακό Μοντέλο
Εδάφους ΕΕ), εφαρμόζει φίλτρα ποιότητας ανά είδος, υπολογίζει εκτιμήσεις
εύρους κατανομής και παράγει τα πινακοποιημένα και χαρτογραφικά αποτελέσματα
που απαιτούνται για την αξιολόγηση κατά το Άρθρο 17.

## Περιεχόμενα

* [Ροή εργασίας](#ροή-εργασίας)
* [Δεδομένα εισόδου](#δεδομένα-εισόδου)
* [Αποτελέσματα](#αποτελέσματα)
* [Κατάλογος ειδών](#κατάλογος-ειδών)
* [Λογισμικό](#λογισμικό)
* [Άδεια χρήσης](#άδεια-χρήσης)

---

## Ροή εργασίας

### 1. Εγκατάσταση

Κλωνοποίηση του αποθετηρίου και επαναφορά του περιβάλλοντος R:

```bash
git clone https://github.com/savvas-paragkamian/necca_epopteia.git
cd necca_epopteia
```

```r
renv::restore()
```

Οι εκδόσεις όλων των πακέτων καταγράφονται στο `renv.lock`. Το `.Rprofile`
ενεργοποιεί αυτόματα το `renv` κατά την εκκίνηση.

#### Container (συνιστάται για πλήρη αναπαραγωγιμότητα)

```bash
podman build -t myproj-rgeo:4.5.3 -f .devcontainer/Containerfile .

podman run --rm -it \
    --userns=keep-id \
    -v "$PWD":/workspaces/project:Z \
    myproj-rgeo:4.5.3
```

Αύξηση πόρων εάν χρειάζεται:

```bash
podman machine stop
podman machine set --memory 12288 --cpus 5
podman machine start
```

### 2. Διαμόρφωση

Όλες οι διαδρομές αρχείων εισόδου και εξόδου δηλώνονται στο `config/params.yml`.
Δεν υπάρχουν σκληρά κωδικοποιημένες διαδρομές στον πηγαίο κώδικα R.

### 3. Εκτέλεση αγωγού

Ο αγωγός υλοποιείται με το πακέτο `{targets}`. Η πλήρης αλυσίδα
Εξαγωγή → Μετασχηματισμός → Φόρτωση εκτελείται με μία εντολή:

```r
targets::tar_make()
```

Το `targets` αποθηκεύει προσωρινά κάθε ενδιάμεσο αποτέλεσμα. Μόνο τα
ξεπερασμένα βήματα επανεκτελούνται όταν αλλάξουν τα δεδομένα εισόδου.
Άλλες χρήσιμες εντολές:

```r
targets::tar_visnetwork()           # οπτικοποίηση γράφου εξαρτήσεων
targets::tar_read(species_range)    # επιθεώρηση αποτελέσματος βήματος
targets::tar_outdated()             # λίστα βημάτων που χρειάζονται επανεκτέλεση
targets::tar_make(par_type = "future")  # παράλληλη εκτέλεση
```

### 4. Δομή αγωγού

Ο αγωγός αποτελείται από 38 στόχους σε τρεις φάσεις.

```
┌─ Εξαγωγή (Extract) ────────────────────────────────────────────────────────┐
│  a17_config · spatial_layers · species_taxonomy                            │
│  10 × read_*_occurrences() (εκτελούνται ανεξάρτητα)                       │
└────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─ Μετασχηματισμός (Transform) ──────────────────────────────────────────────┐
│  national_report_distribution_grid                                         │
│  species_occurrences_invertebrates                                         │
│    → species_samples_art17 → species_samples_eea                           │
│    → species_samples_presence_minimum                                      │
│    → species_samples_presence_elevation                                    │
│    → species_samples_presence_border                                       │
│    → species_samples_presence_dist_flags                                   │
│    → species_samples_presence_pop                                          │
│    → species_samples_presence_final / _final_private                       │
└────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─ Φόρτωση (Load) ───────────────────────────────────────────────────────────┐
│  Χάρτες:  map_borders_quality · map_natura_base · map_art17_overview       │
│           maps_species_occurrences · maps_species_range                    │
│           maps_species_distribution · file_species_range_shp               │
│  TSV:     file_occurrences_invertebrates · file_art17_all                  │
│           file_presence_final · file_distributions_final                   │
│           file_populations_final                                           │
└────────────────────────────────────────────────────────────────────────────┘
```

![Γράφος εξαρτήσεων αγωγού](targets_pipeline_a5_dense.png)

Οι συναρτήσεις R που υλοποιούν κάθε φάση βρίσκονται στον κατάλογο `R/`:

| Αρχείο | Φάση | Ρόλος |
|--------|------|-------|
| `extract_occurrences.R` | Εξαγωγή | Μία συνάρτηση ανάγνωσης ανά πηγή δεδομένων |
| `extract_spatial.R` | Εξαγωγή | Φόρτωση χωρικών επιπέδων αναφοράς |
| `helper_functions.R` | Εξαγωγή/Μετασχηματισμός | Βοηθητικές χωρικές συναρτήσεις (εκτίμηση εύρους, εξαγωγή raster) |
| `transform.R` | Μετασχηματισμός | Εμπλουτισμός, φιλτράρισμα και ανάθεση σημαιών |
| `load_maps.R` | Φόρτωση | Παραγωγή χαρτών ανά είδος και συνολικών χαρτών |
| `load_official_outputs.R` | Φόρτωση | Εγγραφή αρχείων TSV για επίσημη αναφορά Άρθρου 17 |
| `qc.R` | Ποιοτικός έλεγχος | Ποιοτικός έλεγχος (υπό ανάπτυξη) |

---

## Δεδομένα εισόδου

### Δεδομένα παρουσίας ειδών (`data/raw/`)

| Πηγή | Μορφή | Περιεχόμενο |
|------|-------|-------------|
| GBIF | TSV | Καταγραφές, φιλτραρισμένες στην Ελλάδα, αβεβαιότητα συντεταγμένων < 1 000 m |
| E1X ΜΔΠΠ 2014–2024 | XLSX | Δείγματα πεδίου από παρακολούθηση |
| E1X ΒΔ δειγματοληψία | XLSX | Δεδομένα δομημένης δειγματοληψίας |
| E1X ΒΔ βιβλιογραφία | XLSX | Βιβλιογραφικές καταγραφές |
| E2X ΒΔ | TSV | Βάση δεδομένων παρουσίας |
| E2X αναφ. — *Unio crassus* σύμπλεγμα | XLSX | Συμπληρωματικές καταγραφές γλυκόνερων μαλακίων |
| E2X αναφ. — *Stenobothrus eurasius* | TSV | Συμπληρωματικές καταγραφές ακρίδων |
| Ιδιωτικές καταγραφές | CSV | Καταγραφές πεδίου εταίρων έργου |
| Κόκκινος Κατάλογος ΟΦΥΠΕΚΑ | GeoPackage | Σημειακά δεδομένα Κόκκινου Καταλόγου ΟΦΥΠΕΚΑ |
| Εθνική Αναφορά 2013–2018 | SHP | Πλέγμα κατανομής προηγούμενης αναφοράς Άρθρου 17 |
| Σχέδιο Δράσης *Parnassius apollo* 2019 | SHP | Πλέγμα κατανομής για το συγκεκριμένο είδος |

Όλες οι πηγές καταγραφών τυποποιούνται στην ορολογία Darwin Core
(`submittedName`, `decimalLatitude`, `decimalLongitude`, `datasetName`,
`recordNumber`, `collectionCode`, `basisOfRecord`, `individualCount`).

### Χωρικά επίπεδα αναφοράς (`data/spatial/`)

| Επίπεδο | Μορφή | Χρήση |
|---------|-------|-------|
| Διοικητικά όρια Ελλάδας (GADM 4.1) | SHP | Φιλτράρισμα σημείων, υπολογισμός απόστασης από σύνορα |
| Πλέγμα αναφοράς EEA 1 km | SHP | Αξιολόγηση πληθυσμού |
| Πλέγμα αναφοράς EEA 10 km | SHP | Αναφορά κατανομής Άρθρου 17 |
| Natura 2000 (v32, 2021) | SHP | Ανάθεση θέσεων σε περιοχές, χάρτες φόντου |
| Ψηφιακό Μοντέλο Εδάφους ΕΕ (ETRS89-LAEA, 3035) | GeoTIFF | Εξαγωγή υψομέτρου, φίλτρα υψομέτρου ανά είδος |

---

## Αποτελέσματα

Όλες οι διαδρομές εξόδου δηλώνονται στο `config/params.yml`.

| Αρχείο | Τοποθεσία | Περιεχόμενο |
|--------|-----------|-------------|
| `species_occurrences_invertebrates.tsv` | `data/derived/` | Όλες οι καταγραφές παρουσίας ενοποιημένες |
| `species_samples_art17_all.tsv` | `data/derived/` | Καταγραφές φιλτραρισμένες στα είδη Παραρτήματος II |
| `species_samples_presence_final.tsv` | `results/` | Τελικό σύνολο παρουσίας (χωρίς ιδιωτικές καταγραφές) |
| `distributions_presence_final.tsv` | `results/` | Κατανομή ανά είδος και κελί 10 km (παρατηρήσεις + εύρος) |
| `populations_presence_final.tsv` | `results/` | Πληθυσμοί ανά είδος και κελί 1 km |
| `species_range/species_range.shp` | `results/` | Πολύγωνα υπολογισμένου εύρους για όλα τα είδη |
| `maps/map_natura.png` | `results/maps/` | Συνολικός χάρτης Natura 2000 |
| `maps/map_art17_invertebrates_natura.png` | `results/maps/` | Όλες οι καταγραφές Άρθρου 17 σε Natura 2000 |
| `maps/map_occurrences_borders_filtering.png` | `results/maps/` | Χάρτης ποιοτικού ελέγχου: σημεία ως προς τα σύνορα |
| `maps/species_maps/map_*_occurrences.png` | `results/maps/` | Χάρτης καταγραφών ανά είδος |
| `maps/species_maps/map_*_range.png` | `results/maps/` | Χάρτης εύρους κατανομής ανά είδος |
| `maps/species_maps/map_*_distribution.png` | `results/maps/` | Χάρτης κατανομής ανά είδος (πυκνότητα δειγμάτων) |

---

## Κατάλογος ειδών

39 ασπόνδυλα taxa του Παραρτήματος II που αξιολογούνται στον αγωγό:

*Apatura metis* · *Astacus astacus* · *Austropotamobius torrentium* ·
*Bolbelasmus unicornis* · *Buprestis splendens* · *Catopta thrips* ·
*Cerambyx cerdo* · *Coenagrion ornatum* · *Cordulegaster heros* ·
*Dioszeghyana schmidtii* · *Eriogaster catax* · *Euphydryas aurinia* ·
*Euplagia quadripunctaria* · *Hirudo verbana* · *Hyles hippophaes* ·
*Lindenia tetraphylla* · *Lucanus cervus* · *Lycaena dispar* ·
*Maculinea arion* · *Morimus asper funereus* · *Ophiogomphus cecilia* ·
*Osmoderma eremita* Complex · *Papilio alexanor* ·
*Paracaloptenus caloptenoides* · *Parnassius apollo* · *Parnassius mnemosyne* ·
*Polyommatus eroides* · *Probaticus subrugosus* · *Proserpinus proserpina* ·
*Pseudophilotes bavius* · *Rhysodes sulcatus* · *Rosalia alpina* ·
*Stenobothrus eurasius* · *Stylurus flavipes* · *Unio crassus* ·
*Unio elongatulus* · *Vertigo angustior* · *Vertigo moulinsiana* ·
*Zerynthia polyxena*

---

## Λογισμικό

R 4.5.3 διαχειριζόμενο μέσω `renv`. Κύρια πακέτα:

| Πακέτο | Έκδοση | Ρόλος |
|--------|--------|-------|
| targets | 1.12.0 | Ενορχήστρωση αγωγού |
| tarchetypes | 0.14.1 | Βοηθητικά εργαλεία αγωγού |
| sf | 1.0-19 | Διανυσματικά χωρικά δεδομένα |
| terra | 1.8-15 | Raster χωρικά δεδομένα |
| dplyr / tidyr / purrr | tidyverse | Επεξεργασία δεδομένων |
| ggplot2 | 3.5.1 | Χαρτογράφηση |
| readxl | 1.4.3 | Ανάγνωση Excel |
| rgbif | 3.8.1 | Πρόσβαση σε δεδομένα GBIF |
| taxize | 0.10.0 | Ταξινομική εναρμόνιση |
| yaml | — | Ανάλυση αρχείου διαμόρφωσης |

---

## Άδεια χρήσης
MIT License
Παραγκαμιάν, Σ., Τζωρτζακάκη, Ό., Γουδέλη, Γ., Βούρκα, Α. και Κασσάρα, Χ., ΟΦΥΠΕΚΑ, Ελλάδα
