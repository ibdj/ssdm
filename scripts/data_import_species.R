#### packages ##################################################################

library(tidyverse)
library(vegan)
library(janitor)
library(terra)
library(sf)

#### functions #################################################################

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

generate_dataframe <- function(number, data = survey_0) {
  taxon_col  <- paste0("taxon_", number)
  height_col <- paste0("taxon_", number, "_height")
  bb_col     <- paste0("taxon_", number, "_bb")
  
  # everything that is NOT a per-taxon column = the shared plot metadata
  taxon_pattern <- "^taxon_\\d+(_height|_braun_blanquet)?$"
  meta_cols <- setdiff(names(data), grep(taxon_pattern, names(data), value = TRUE))
  
  data |>
    dplyr::select(dplyr::all_of(c(meta_cols, taxon_col, height_col, bb_col))) |>
    dplyr::mutate(position = paste0("taxon_", number)) |>
    dplyr::rename(taxon  = dplyr::all_of(taxon_col),
                  height = dplyr::all_of(height_col),
                  bb     = dplyr::all_of(bb_col))
}

#### importing survey123 file ##################################################

survey_0 <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/survey_0.csv", col_types = cols(Date = col_datetime(format = "%m/%d/%Y %H.%M"))) |> 
  dplyr::select(where(~ !all(is.na(.))))

survey_0 <- survey_0 |> 
  clean_names() |>  
  mutate(rowid = row_number()) |> 
  dplyr::rename_with(~ gsub("braun_blanquet", "bb", .x)) 

# find the highest taxon slot present in the data
n_taxa <- names(survey_0) |>
  stringr::str_extract("(?<=^taxon_)\\d+") |>   # pull the number after "taxon_"
  as.integer() |>
  max(na.rm = TRUE)



raw_df <- generate_dataframe(survey_0)
  
names(survey_0)
str(survey_0)
summary(survey_0)

#### species matrix ############################################################

# Step 1: pivot just the species names to long
taxon_names <- raw_df_cover |>
  dplyr::select(plot_name, matches("^taxon_[0-9]+$")) |>
  pivot_longer(-plot_name, names_to = "slot", values_to = "species_name")

# Step 2: pivot just the bb values to long
taxon_bb <- raw_df_cover |>
  dplyr::select(plot_name, matches("^taxon_[0-9]+_bb$")) |>
  pivot_longer(-plot_name, names_to = "slot", values_to = "cover") |>
  mutate(slot = str_remove(slot, "_bb$"))

# Step 3: pivot just the height values to long
taxon_height <- raw_df_cover |>
  dplyr::select(plot_name, matches("^taxon_[0-9]+_height$")) |>
  pivot_longer(-plot_name, names_to = "slot", values_to = "height") |>
  mutate(slot = str_remove(slot, "_height$"))

# Step 4: join all three together
species_long <- taxon_names |>
  left_join(taxon_bb, by = c("plot_name", "slot")) |>
  left_join(taxon_height, by = c("plot_name", "slot")) |>
  filter(!is.na(species_name) & species_name != "") |> 
  dplyr::select(-slot)

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
  distinct(species_name) |> 
  arrange(species_name) |> 
  print(n = Inf)

species_long |> 
  dplyr::count(plot_name, species_name) |> 
  dplyr::filter(n > 1)

species_matrix <- species_long |>
  dplyr::select(plot_name, species_name, cover) |>
  pivot_wider(
    names_from = species_name,
    values_from = cover,
    values_fill = 0
  )

sp_cols <- species_matrix |> dplyr::select(-plot_name)
