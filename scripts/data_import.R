#### packages ####

library(tidyverse)
library(janitor)

#### loading the data ####

tms <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/tms_pivot.rds") |> 
  clean_names() |> 
  group_by(plot) |> 
  reframe(temp_mean_tms = mean(temp),
          mean_soilmoisture_tms = mean(moisture_data_raw))

summary(tms)

samples_qgis <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/samples_qgis.csv") |> 
  select(plot, X,Y,elevation, ndvi, ndwi) |> 
  left_join(tms, by = "plot")

names(samples_qgis)

samples <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/samples.csv")

names(samples)

