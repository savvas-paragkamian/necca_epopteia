# Analysis of the invertebrates of Greece


## IUCN Greece 

All the data available from iunc Search on 2024-10-30 at 08:47:29

```
gawk -F"," '(NR>1){a[$3]++}END{for (i in a){print i "\t" a[i]}}' points_data.csv | sort -n -k 3
```

614 species with points data with 753110 occurrences.

## NECCA database schema

It is under development from AUTH.


## Invertebrate species defined in Natura2000 

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

## Natura2000

The MS access database has the schema of EEA Natura2000 in access format.

To retrieve each table run the following 

```
mdb-export -d "\t" Natura2000DB_V32.mdb other_species > Natura2000DB_V32_other_species.tsv
```

The [EEA Natura2000 data](https://www.eea.europa.eu/data-and-maps/data/natura-14/natura-2000-tabular-data-12-tables/#SPECIES).

## Data retrieval

### GBIF database 

Here we retrieve all the occurrences of the species from GBIF.



## Spatial analysis

### Favourable Reference Values - FRVs

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


### DEM

* digital elevation model over Europe (EU-DEM) [eurostat DEM](https://ec.europa.eu/eurostat/web/gisco/geodata/digital-elevation-model/eu-dem#DD)
* slope, aspect


