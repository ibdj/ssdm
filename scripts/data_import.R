#### packages ####

library(tidyverse)
library(vegan)
library(janitor)

#### functions ####
bb_to_cover <- function(x) {
  dplyr::case_when(
    startsWith(x, "5 (") ~ 87.5,
    startsWith(x, "4 (") ~ 62.5,
    startsWith(x, "3 (") ~ 37.5,
    startsWith(x, "2 (") ~ 15.0,
    startsWith(x, "1 (") ~ 2.5,
    startsWith(x, "+ (") ~ 1.0,
    startsWith(x, "r (") ~ 0.5,
    startsWith(x, "i (") ~ 0.1,
    x == "0" ~ 0.0,
    TRUE ~ NA_real_
  )
}

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

df_cover <- df_raw |> 
 mutate(across(ends_with("_bb"), bb_to_cover)) |> 
  mutate(total_cover = rowSums(across(ends_with("_bb")), na.rm = TRUE),
         richness = rowSums(across(ends_with("_bb")) > 0, na.rm = TRUE))
  

abiotic_plot <- df_cover |> 
  dplyr::select(plot_name , veg_height_ave, bare_ground_bb, x, y, total_cover, richness, shannon)
  
summary(abiotic_plot)


species_cols <- df_raw |> 
  select(ends_with("_bb")) |> 
  names() |> 
  str_remove("_bb$") |> 
  keep(~ .x %in% names(df_raw))  # only keep if taxon column actually exists

# Then run the mismatch check as before
mismatches <- map_dfr(species_cols, function(sp) {
  df_raw |> 
    filter(!is.na(.data[[sp]]) & .data[[sp]] != "" & 
             is.na(.data[[paste0(sp, "_bb")]])) |> 
    mutate(species = sp, 
           taxon_value = .data[[sp]],
           bb_value = .data[[paste0(sp, "_bb")]]) |> 
    select(rowid, species, taxon_value, bb_value)
})

mismatches

n_taxa <- df_raw |>
  select(matches("^taxon_\\d+$")) |>
  ncol()

df <- map(seq_len(n_taxa), \(number) {
  taxon_col  <- sym(paste0("taxon_", number))
  height_col <- sym(paste0("taxon_", number, "_height"))
  bb_col     <- sym(paste0("taxon_", number, "_bb"))
  
  df_raw |>
    select(1:31, !!taxon_col, !!height_col, !!bb_col, 74:78) |>
    mutate(position = paste0("taxon_", number)) |>
    rename(taxon  = !!taxon_col,
           height = !!height_col,
           bb     = !!bb_col)
}) |>
  list_rbind()

df <- df |> 
  mutate(taxon = as.factor(taxon))

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
