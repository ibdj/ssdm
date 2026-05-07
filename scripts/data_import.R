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
  mutate(total_cover = rowSums(across(ends_with("_bb")), na.rm = TRUE))

summary(df_cover)
  
#### final abiotic df###########################################################

abiotic_plot <- df_cover |> 
  dplyr::select(plot_name , veg_height_ave, bare_ground_bb, x, y, total_cover, richness, shannon)
  
summary(abiotic_plot)


#### species matrix ############################################################

# Step 1: pivot just the species names to long
taxon_names <- df_cover |>
  select(plot_name, matches("^taxon_[0-9]+$")) |>
  pivot_longer(-plot_name, names_to = "slot", values_to = "species_name")

# Step 2: pivot just the bb values to long
taxon_bb <- df_cover |>
  select(plot_name, matches("^taxon_[0-9]+_bb$")) |>
  pivot_longer(-plot_name, names_to = "slot", values_to = "cover") |>
  mutate(slot = str_remove(slot, "_bb$"))

# Step 3: pivot just the height values to long
taxon_height <- df_cover |>
  select(plot_name, matches("^taxon_[0-9]+_height$")) |>
  pivot_longer(-plot_name, names_to = "slot", values_to = "height") |>
  mutate(slot = str_remove(slot, "_height$"))

# Step 4: join all three together
species_long <- taxon_names |>
  left_join(taxon_bb, by = c("plot_name", "slot")) |>
  left_join(taxon_height, by = c("plot_name", "slot")) |>
  filter(!is.na(species_name) & species_name != "") |> 
  select(-slot)

species_long |> 
  distinct(species_name) |> 
  arrange(species_name) |> 
  print(n = Inf)

species_long <- species_long |>
  mutate(
    species_name = str_trim(species_name),           # remove whitespace
    species_name = str_remove(species_name, "_+$"),  # remove trailing underscores
    species_name = case_when(
      species_name == "Scirpis caespitosus" ~ "Scirpus caespitosus",
      TRUE ~ species_name
    )
  ) |> 
  group_by(plot_name, species_name) |>
  slice_max(cover, n = 1, with_ties = FALSE) |>
  ungroup()

species_long |> 
  count(plot_name, species_name) |> 
  filter(n > 1)

species_matrix <- species_long |>
  select(plot_name, species_name, cover) |>
  pivot_wider(
    names_from = species_name,
    values_from = cover,
    values_fill = 0
  )

sp_cols <- species_matrix |> select(-plot_name)

abiotic_plot <- abiotic_plot |>
  left_join(species_matrix |> select(plot_name), by = "plot_name") |>
  mutate(
    richness = rowSums(sp_cols > 0),
    shannon = vegan::diversity(sp_cols, index = "shannon")
  )

#### species frequency #########################################################

species_frequency <- species_matrix |>
  select(-plot_name) |>
  summarise(across(everything(), ~ sum(. > 0))) |>
  pivot_longer(everything(), 
               names_to = "species", 
               values_to = "n_plots") |>
  arrange(n_plots)

print(species_frequency, n = Inf)

#### other stuff ###############################################################
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
