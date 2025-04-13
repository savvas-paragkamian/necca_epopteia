# Analysis of the invertebrates of Greece

This repository contains the scripts for the analysis and visualisation
of the monitoring of invertebrates as part of the Natura2000 Species 
directive (Art. 17 under the directive 43/92/ECC). 

This work is funded by the Natural Environment & Climate Change Agency during 
the period 10-2024 until 07-2025.

## Contents

* [Scripts](#scripts)
* [Species enrichment](#species-enrichment)
* [Spatial analysis](#spatial-analysis)
* [Favourable Reference Values](#favourable-reference-values)
* [Natura2000 targets](#natura2000-targets)
* [Data validation and quality control](#data-validation-and-quality-control)
* [NECCA database](#necca-database)
* [Software](#software)
* [Licence](#licence)

## Scripts

The scripts in this repository are (with the order of execution):
1. species_enrichment.R, species info gathering and homogenisation
2. spatial_analysis.R, spatial data trimming and final enrichment of species occurrences
3. frvs_invertebrates_gr.R
4. N2000_targets_invertebrates_gr.R

The script *data_validation_qc.R* is for the data validation of 
the data from the current samplings in order to find errors prior to
submission.

The script *necca_db_standards.R* aims to construct all the data 
compliled from this work in a compatible form with the NECCA database.

## Species enrichment

### Table with invertebrates information

| Name                          | Format            | Content               | Source                              |
|-------------------------------|------------------|------------------------|-------------------------------------|
| Species Taxonomy              | Tabular (xlsx)   | Species info           | Manual curation - EIONET            |
| Directive value               | Tabular (xlsx)   | Species info           | ENVENCO                             |
| Current value                 | Tabular (xlsx)   | Species info           | Παραδοτέο Ανάδοχου                  |
| IUCN Red list                 | Tabular (tsv)    | Species info           | IUCN                                |
| IUCN Red list                 | Tabular (csv)    | Species occurrences    | IUCN                                |
| NECCA Redlist                 | Tabular (xlsx)   | Species info           | NECCA db                            |
| NECCA Redlist                 | Tabular (xlsx)   | Species occurrences    | NECCA db                            |
| Art 17 1.3.1                  | Tabular (xlsx)   | Species info           | EIONET                              |
| Natura2000_v32                | Access           | Species info           | EIONET                              |
| Gbif                          | Tabular (tsv)    | Species occurrences    | Gbif                                |
| ΕΠΟΠΤΕΙΑ Ι                    | Tabular (xlsx)   | Species occurrences    | ENVENCO                             |
| ΕΠΟΠΤΕΙΑ Ι                    | Tabular (xlsx)   | Species info           | ENVENCO                             |
| ΟΦΥΠΕΚΑ έργα Μονάδων          | Tabular (xlsx)   | Species occurrences    | NECCA db                            |
| Edaphobase 2.0                | Tabular (csv)    | Species occurrences    | Edaphobase 2.0                      |
| National Report 2013-2018     | Shp              | Species info           | ENVENCO                             |
| Molluscabase                  | ?                | Species occurrences    | [molluscabase.org](https://www.molluscabase.org/) |
| Worms Traits                  | Tabular (tsv)    | Species info           |                                     |
| Genetic diversity             |                  |                        |                                     |
| Species threats?              |                  |                        |                                     |

### Invertebrate species of Annex II of article 17 for the Natura2000 

Apatura metis
Astacus astacus
Austropotamobius torrentium
Bolbelasmus unicornis
Buprestis splendens
Catopta thrips
Cerambyx cerdo
Coenagrion ornatum
Cordulegaster heros
Dioszeghyana schmidtii
Eriogaster catax
Euphydryas aurinia
Euplagia quadripunctaria
Hyles hippophaes
Lindenia tetraphylla
Lucanus cervus
Lycaena dispar
Maculinea arion
Ophiogomphus cecilia
Papilio alexanor
Paracaloptenus caloptenoides
Parnassius apollo
Parnassius mnemosyne
Polyommatus eroides
Probaticus subrugosus
Proserpinus proserpina
Pseudophilotes bavius
Rhysodes sulcatus
Rosalia alpina
Stenobothrus eurasius
Stylurus flavipes
Unio crassus
Unio elongatulus
Vertigo angustior
Vertigo moulinsiana
Zerynthia polyxena
Morimus asper funereus
Osmoderma eremita Complex
Hirudo verbana


Endemic species and critically endangered are also studied in order to include
them in future analyses.



### IUCN Greece 

All the data available from iunc Search on 2024-10-30 at 08:47:29

```
gawk -F"," '(NR>1){a[$3]++}END{for (i in a){print i "\t" a[i]}}' points_data.csv | sort -n -k 3
```

614 species with points data with 753110 occurrences.


### GBIF database 

Here we retrieve all the occurrences of the species from GBIF.


## Spatial analysis

### Table with spatial data

| Name                             | Format| Content                          | Source                        |
|----------------------------------|-------|-------------------------------------|-----------------------------|
| Natura 2000 Habitat types        | shp   | habitat                             | Κτηματολόγιο (ΕΚΧΑ)        |
| Natura 2000 vegetation types     | shp   | vegetation                          | Διεύθυνση δασών            |
| Corine Land Cover                | shp   | Land use                            | EEA                         |
| Natura2000                       | shp   | Natura2000 polygons                 | EIONET                      |
| EAA grid                         | shp   | 1km², 10km²                          | EEA                         |
| Digital Elevation model         | tiff  |                                     | EEA                         |
| Slope                            |       |                                     |                             |
| Exposure                         |       |                                     |                             |
| Springs                          |       | ?                                   | ?                           |
| WorldClim                        | tiff  |                                     | WorldClim 2.0              |
| GlobCover                        | ?     |                                     | ?                           |
| HildaPlus                        | tiff  | Land cover change from 1960         |                             |
| Gadm                             | shp   | Administrative units                | Gadm                        |
| Υδρογραφικό δίκτυο Ελλάδα       | shp   |                                     | ?                           |
| Natura2000 threats               | ?     |                                     |                             |
| Elstat                           |       |                                     |                             |

### Natura2000

The MS access database has the schema of EEA Natura2000 in access format.

To retrieve each table run the following 

```
mdb-export -d "\t" Natura2000DB_V32.mdb other_species > Natura2000DB_V32_other_species.tsv
```

The [EEA Natura2000 data](https://www.eea.europa.eu/data-and-maps/data/natura-14/natura-2000-tabular-data-12-tables/#SPECIES).

### DEM

* digital elevation model over Europe (EU-DEM) [eurostat DEM](https://ec.europa.eu/eurostat/web/gisco/geodata/digital-elevation-model/eu-dem#DD)
* slope, aspect

### Slope

## Favourable Reference Values

The paper: *Defining applying the concept of Favourable Reference Values for species and habitats*
defines the FRVs.

### MaxEnt

The following variables are used mostly: 
* Altitude
* Slope
* Aspect
* Annual Mean Temperature
* Temperature Seasonality
* Mean Temp of Wettest Quarter
* Annual Precipitation
* Precipitation Seasonality
* Precipitation of Warmest Quarter

## Natura2000 targets

## Data validation and quality control

## NECCA database

It is under development from AUTH. There is a script under development for the 
standardisation of data in order to be uploaded to NECCA database.

## Software

R version 4.4.3 (2025-02-28) with the packages:
* tidyr_1.3.1
* tibble_3.2.1
* ggplot2_3.5.1 
* readr_2.1.5
* purrr_1.0.2 
* dplyr_1.1.4
* readxl_1.4.3
* rgbif_3.8.1
* vegan_2.6-10
* taxize_0.10.0
* terra_1.8-15
* sf_1.0-19

GNU Awk 5.3.1

## Licence

GNU GPLv3 license (for 3rd party scripts separate licenses apply).
