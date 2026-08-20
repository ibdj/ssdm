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

eco_veg_growth_forms <- read_excel("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/01 Vegetation changes Kobbefjord/data/nero_analysis/data/eco_veg_growth_forms.xlsx") |> 
  mutate(taxon = species, func_type = as.factor(func_type))
names(eco_veg_growth_forms)

survey_0 <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/survey_0.csv", col_types = cols(Date = col_datetime(format = "%m/%d/%Y %H.%M"))) |> 
  dplyr::select(where(~ !all(is.na(.)))) |> 
  clean_names() |>  
  mutate(rowid = row_number()) |> 
  dplyr::rename_with(~ gsub("braun_blanquet", "bb", .x)) |> 
  mutate(plot_name = toupper(plot_name)) |> 
  filter(grepl('MP', plot_name))
  

# find the highest taxon slot present in the data
n_taxa <- names(survey_0) |>
  stringr::str_extract("(?<=^taxon_)\\d+") |>   # pull the number after "taxon_"
  as.integer() |>
  max(na.rm = TRUE)

names(survey_0)
str(survey_0)
summary(survey_0)

survey_0 <- survey_0 |> 
  dplyr::select(
    plot_name,
    vegetation_height_n,
    vegetation_height_e,
    vegetation_height_s,
    vegetation_height_w,
    soil_moisture_n,
    soil_moisture_e,
    soil_moisture_s,
    soil_moisture_w,
    soil_temp_n,
    soil_temp_e,
    soil_temp_s,
    soil_temp_w,
    topographic_complexity_cm,
    bryophyte_bb_49,
    lichen_bb,
    other_notes,
    vegetation_type,
    bare_ground_bb,
    other_vegetation_type,
    tms,
    x,
    y,
    rowid,
    taxon_1,
    taxon_1_height,
    taxon_1_bb,
    taxon_2,
    taxon_2_height,
    taxon_2_bb,
    taxon_3,
    taxon_3_height,
    taxon_3_bb,
    taxon_4,
    taxon_4_height,
    taxon_4_bb,
    taxon_5,
    taxon_5_height,
    taxon_5_bb,
    taxon_6,
    taxon_6_height,
    taxon_6_bb,
    taxon_7,
    taxon_7_height,
    taxon_7_bb,
    taxon_8,
    taxon_8_height,
    taxon_8_bb,
    taxon_9,
    taxon_9_height,
    taxon_9_bb,
    taxon_10,
    taxon_10_height,
    taxon_10_bb,
    taxon_11,
    taxon_11_height,
    taxon_11_bb,
    taxon_12,
    taxon_12_height,
    taxon_12_bb,
    taxon_13,
    taxon_13_height,
    taxon_13_bb,
    taxon_14,
    taxon_14_height,
    taxon_14_bb,
    taxon_15,
    taxon_15_height,
    taxon_15_bb,
    taxon_16,
    taxon_16_height,
    taxon_16_bb
  )

survey_0_ren <- survey_0 |>
  dplyr::rename_with(
    ~ paste0(.x, "_name"),
    .cols = dplyr::matches("^taxon_\\d+$")   # $ = ONLY bare taxon_N, not _height/_bb
  )

species_long <- survey_0_ren |>
  tidyr::pivot_longer(
    cols = dplyr::matches("^taxon_\\d+"),
    names_to = c("slot", ".value"),
    names_pattern = "taxon_(\\d+)_(name|height|bb)"
  ) |>
  dplyr::filter(!is.na(name)) |>
  dplyr::rename(taxon = name) |> 
  mutate(
    taxon = str_trim(taxon),
    taxon = str_remove(taxon, "_+$"),
    taxon = case_when(
      taxon == "Scirpis caespitosus" ~ "Scirpus caespitosus",
      TRUE ~ taxon
    )
  )

species_long <- species_long |>
  dplyr::mutate(cover = bb_to_cover(bb)) |> 
  group_by(plot_name, taxon) |>
  slice_max(cover, n = 1, with_ties = FALSE) |>
  ungroup() |> 
  left_join(eco_veg_growth_forms, by = "taxon")

taxon_check <- species_long |> 
  filter(is.na(ecoveg_sgfc)) |> 
  distinct(taxon)
#### species matrix ############################################################

species_only_long <- species_long |> 
  dplyr::select(plot_name, taxon, cover) |> 
  distinct()

species_matrix <- species_long |>
  dplyr::select(plot_name, taxon, cover) |>
  tidyr::pivot_wider(names_from = taxon, values_from = cover, values_fill = 0) |>
  dplyr::right_join(dplyr::distinct(survey_0, plot_name), by = "plot_name") |>
  dplyr::mutate(dplyr::across(-plot_name, ~ tidyr::replace_na(.x, 0)))

species_matrix_cover <- species_matrix |>
  tibble::column_to_rownames("plot_name")

species_matrix_pa <- (species_matrix_cover > 0) * 1

saveRDS(species_matrix_out, "data/species_matrix_cover.rds")
saveRDS(species_matrix_pa, "data/species_matrix_pa.rds")


