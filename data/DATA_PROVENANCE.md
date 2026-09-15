# Data provenance and licence notes

This repository contains selected data needed for the post-hoc examples. Data files are not relicensed by this repository; they retain the licences, terms, and attribution requirements of the original providers.

## `data/no2_de/no2_2018.gpkg`

Derived from German Umweltbundesamt (UBA) station metadata and annual mean NO2 observations for 2018. The UBA Luftdaten portal provides current and historical air-quality measurements and an API for automated data access. The UBA site states that federal geodata and geodata services are subject to the German GeoNutzV terms. Cite the UBA data sources in the manuscript and check the current UBA terms before redistribution.

## `data/no2_eu/NO2_2010_prox.csv` and coordinate-augmented files

The station table originates from the data accompanying Vizcaino and Lavalle (2018), *Development of European NO2 Land Use Regression Model for present and future exposure assessment: Implications for policy analysis*, Mendeley Data, version 1, DOI [`10.17632/kkss4z9yrs.1`](https://doi.org/10.17632/kkss4z9yrs.1). The Mendeley Data record states that the dataset is licensed under the Creative Commons Attribution 4.0 International licence (CC BY 4.0). The coordinate-augmented files are modified derivatives that add station coordinates recovered from historical AirBase/European Environment Agency metadata; attribution to Vizcaino and Lavalle and an indication of these changes must be retained.

## `data/airbase_metadata/`

Historical AirBase/EEA station metadata used to recover coordinates for the NO2--EU example. The EEA AirBase export page documents metadata exports and links to the Pan-European metadata CSV. Check current EEA reuse terms before redistributing modified versions.

## Built-in and reported-data examples

The Meuse examples use `sp::meuse` and `sp::meuse.grid`. Biomass, BTS, SWE, S1--SM, and SMAP--CRNS use numerical values reported in the cited papers; no raw source data from those papers are redistributed here.
