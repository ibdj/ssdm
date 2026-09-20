library(rgbif)
library(CoordinateCleaner)
library(sf)
library(terra)
library(dplyr)

my_species <- read_rds("data/species_frequency.rds") |> 
  select(taxon)

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

# WWF ecoregions: download "official teow" shapefile (e.g. via WWF site)
eco <- st_read("/Users/ibdj/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/wwf_terr_ecos/wwf_terr_ecos.shp") |>
  filter(BIOME %in% c(6, 11)) |>          # 6 = boreal, 11 = tundra
  st_transform("EPSG:6931")

grid <- rast(ext(eco), resolution = 50000, crs = "EPSG:6931")
biome_r <- rasterize(vect(eco), grid, field = "BIOME")

biome_freq <- terra::freq(biome_r)
n_cells_boreal_total <- biome_freq$count[biome_freq$value == 6]
n_cells_tundra_total <- biome_freq$count[biome_freq$value == 11]

# 4b. Assign each occurrence to a cell + biome  -> occ_cells
pts <- occ |>
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) |>
  st_transform("EPSG:6931") |>
  vect()

occ_cells <- occ |>
  mutate(cell  = cells(biome_r, pts)[, "cell"],
         biome = terra::extract(biome_r, pts)[, 2]) |>
  filter(!is.na(biome)) |>                # drops sea, temperate, everything else
  mutate(biome = if_else(biome == 6, "boreal", "tundra"))

# 5. Index: occupied 50 km cells per species per biome
idx <- occ_cells |>
  dplyr::distinct(species, cell, biome) |>
  dplyr::summarise(n = dplyr::n(), .by = c(species, biome)) |>
  tidyr::pivot_wider(names_from = biome, values_from = n, values_fill = 0) |>
  dplyr::mutate(B = boreal / (boreal + tundra), 
                B_abs = (boreal / n_cells_boreal_total) /
                  (boreal / n_cells_boreal_total + tundra / n_cells_tundra_total))

