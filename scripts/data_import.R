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

df_raw <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/samples.csv", col_types = cols(Date = col_datetime(format = "%m/%d/%Y %H.%M"))) |> 
  clean_names() |> 
  mutate(rowid = row_number())


generate_dataframe <- function(number) {
  taxon_col <- sym(paste0("taxon_", number))
  height_col <- sym(paste0("taxon_", number, "_height"))
  bb_col <- sym(paste0("taxon_", number, "_bb"))
  
  df_raw %>%
    select(1:31, !!taxon_col, !!height_col, !!bb_col, 74:78) %>%
    mutate(rowid = row_number(),
           position = paste0("taxon_", number)) %>%
    rename(taxon = !!taxon_col,
           height = !!height_col,
           bb = !!bb_col)
}

taxon_list <- lapply(1:14, generate_dataframe)

pivot <- bind_rows(taxon_list) %>% 
  mutate(veg_mean_height = rowMeans(select(.,veg_height_n,veg_height_s,veg_height_e,veg_height_w)))%>%
  filter(!is.na(taxon))
