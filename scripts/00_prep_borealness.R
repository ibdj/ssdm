library(rgbif)
library(CoordinateCleaner)
library(sf)
library(terra)
library(dplyr)

# 1. Match your names to the GBIF backbone
keys <- name_backbone_checklist(my_species) |>
  filter(matchType %in% c("EXACT", "FUZZY")) |>   # inspect fuzzy ones manually
  pull(usageKey)

# 2. One bulk download (needs GBIF credentials in .Renviron)
d <- occ_download(
  pred_in("taxonKey", keys),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_lt("coordinateUncertaintyInMeters", 10000),
  pred_in("basisOfRecord", c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION")),
  format = "SIMPLE_CSV"
)
occ_download_wait(d)
occ <- occ_download_get(d) |> occ_download_import()

# 3. Clean
occ <- occ |> clean_coordinates(lon = "decimalLongitude", lat = "decimalLatitude",
                                tests = c("centroids", "institutions", "equal", "zeros", "seas")) |>
  filter(.summary)

# 4. Biome raster (once): WWF ecoregions -> boreal (6) / tundra (11) on EPSG:6931 grid
# 5. Index: occupied 50 km cells per species per biome
idx <- occ_cells |>
  distinct(species, cell, biome) |>
  count(species, biome) |>
  tidyr::pivot_wider(names_from = biome, values_from = n, values_fill = 0) |>
  mutate(B = boreal / (boreal + tundra))