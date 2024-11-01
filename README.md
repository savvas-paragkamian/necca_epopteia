# Analysis of the invertebrates of Greece


## IUCN Greece 

All the data available from iunc Search on 2024-10-30 at 08:47:29

```
gawk -F"," '(NR>1){a[$3]++}END{for (i in a){print i "\t" a[i]}}' points_data.csv | sort -n -k 3
```

614 species with points data with 753110 occurrences.



## Favourable Reference Values - FRVs

The paper: *Defining applying the concept of Favourable Reference Values for species and habitats*
defines the FRVs.


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
