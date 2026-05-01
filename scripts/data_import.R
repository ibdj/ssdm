#### packages ####

library(tidyverse)

#### loading the data ####

tms <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/tms_pivot.rds")

names(tms)

tms_means <- tms |> 
  group_by(plot) |> 
  reframe(temp_mean = mean(temp),
         mean_soilmoisture = mean(moisture_data_raw))

samples_qgis <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/samples_qgis.csv")

names(samples_qgis)

samples <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/samples.csv")
