#### packages ##################################################################

library(tidyverse)
library(vegan)
library(janitor)
library(terra)
library(sf)
library(whitebox) # to calculate twi
whitebox::install_whitebox()  # installs the WhiteboxTools binary
library(readxl)
library(gstat)
#devtools::install_github('cjcarlson/embarcadero') #sdm also for interpolation
library(embarcadero)
library(pROC)
library(dbarts)
library(raster)
library(caret)
library(GGally) # to make correlation plots for the soil moisture
library(spatialEco) # for radiation / heat load index
library(car)

#### cover matrix ##############################################################
# Build cover matrix (replaces pa_matrix)
cover_matrix <- species_matrix |>
  dplyr::select(plot_name, all_of(modelable_species)) |>
  left_join(
    mp_abiotic |>
      dplyr::select(plot_name, elevation, slope, hli, ndwi, temp, snowfree),
    by = "plot_name"
  )