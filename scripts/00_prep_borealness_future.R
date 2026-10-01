library(rgbif)
library(CoordinateCleaner)
library(sf)
library(terra)
library(tidyverse)
library(ggplot2)
library(janitor)
library(readxl)
library(tidyterra)

modelable_species <- read_rds("data/modelable_species.rds")

#my_species <- read_rds("data/species_frequency.rds") |> 
#  select(taxon)

my_species <- as_tibble(modelable_species) |> 
  mutate(scientificName = value) |> 
  select(scientificName)

matched_modelable <- purrr::map(my_species$scientificName, \(x) {
  r <- tryCatch(rgbif::name_backbone(name = x), error = \(e) NULL)
  if (is.null(r)) { Sys.sleep(5); r <- rgbif::name_backbone(name = x) }
  r
}, .progress = TRUE) |>
  purrr::list_rbind()

saveRDS(matched_modelable, "data/gbif_backbone_match_modelable.rds")

# 1. Match your names to the GBIF backbone
keys_modelable <- matched_modelable |>
  filter(matchType %in% c("EXACT", "FUZZY")) |>   # inspect fuzzy ones manually
  pull(usageKey)

saveRDS(keys_modelable, "data/keys_modelable.rds")

# 2. One bulk download (needs GBIF credentials in .Renviron)
d_modelable <- occ_download(
  pred_in("taxonKey", keys_modelable),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_lt("coordinateUncertaintyInMeters", 10000),
  pred_in("basisOfRecord", c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION")),
  format = "SIMPLE_CSV"
)

occ_download_wait(d_modelable)
occ_modelable <- occ_download_get(d_modelable) |> occ_download_import()
saveRDS(occ_modelable, "data/occ_raw_modelable.rds")

z <- occ_download_get(d, overwrite = TRUE)   # fresh download of the zip

occ <- data.table::fread(
  cmd = "unzip -p 0009781-260921141020460.zip",
  select = c("species", "decimalLongitude", "decimalLatitude",
             "coordinateUncertaintyInMeters", "countryCode",
             "basisOfRecord", "year", "issue"),
  quote = ""
)

dplyr::n_distinct(occ$species)     # should now be ~460-480, not 74
saveRDS(occ, "data/occ_raw_full.rds")sa

# 3. Clean
occ_clean <- occ |>
  as.data.frame() |>
  cc_val(lon = "decimalLongitude", lat = "decimalLatitude") |>
  cc_equ(lon = "decimalLongitude", lat = "decimalLatitude") |>
  cc_zero(lon = "decimalLongitude", lat = "decimalLatitude") |>
  cc_cen(lon = "decimalLongitude", lat = "decimalLatitude", buffer = 1000) |>
  cc_inst(lon = "decimalLongitude", lat = "decimalLatitude")

saveRDS(occ_clean, "data/occ_clean.rds")

occ_clean <- readRDS("data/occ_clean.rds")


# chekcing the lenght of the species lists ###

length(keys)                              # ~500 expected
dplyr::n_distinct(occ$species)            # the download
dplyr::n_distinct(occ_clean$species)
dplyr::n_distinct(occ_cells$species)
dplyr::n_distinct(idx$species)

# 4. Biome raster (once): WWF ecoregions -> boreal (6) / tundra (11) on EPSG:6931 grid

# WWF ecoregions: download "official teow" shapefile (e.g. via WWF site)

#---- Biome raster from WWF ecoregions (boreal = 6, tundra = 11) ----

eco_raw <- st_read("/Users/ibdj/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/wwf_terr_ecos/wwf_terr_ecos.shp")

eco <- eco_raw |>
  dplyr::filter(BIOME %in% c(6, 11)) |>
  st_make_valid()

# drop Southern Hemisphere polygons (Antarctic/Patagonian "tundra")
eco <- eco[sapply(st_geometry(eco), \(g) st_bbox(g)["ymax"] > 30), ] |>
  st_transform("EPSG:6931")

stopifnot(nrow(eco) > 0)   # guard: fail loudly, not with an empty raster

# 50 km equal-area grid
eco_v   <- vect(eco)
grid    <- rast(ext(eco_v), resolution = 50000, crs = "EPSG:6931")
biome_r <- rasterize(eco_v, grid, field = "BIOME")
writeRaster(biome_r, "data/biome_r.tif")

# checking that all data is save
file.exists(c("data/gbif_backbone_match.rds", "data/occ_raw.rds", "data/occ_clean.rds", "data/biome_r.tif", "data/keys.rds"))

biome_r <- rast("data/biome_r.tif")     # if the tif exists

biome_freq <- terra::freq(biome_r)
n_cells_boreal_total <- biome_freq$count[biome_freq$value == 6]
n_cells_tundra_total <- biome_freq$count[biome_freq$value == 11]

biome_freq <- terra::freq(biome_r)
n_cells_boreal_total <- biome_freq$count[biome_freq$value == 6]
n_cells_tundra_total <- biome_freq$count[biome_freq$value == 11]
n_cells_boreal_total / n_cells_tundra_total   # sanity check ratio

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
                B_std = (boreal / n_cells_boreal_total) /
                  (boreal / n_cells_boreal_total + tundra / n_cells_tundra_total))

idx |>
  tidyr::pivot_longer(c(B, B_std), names_to = "index", values_to = "value") |>
  ggplot(aes(value)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "grey30") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  facet_wrap(~index, labeller = as_labeller(
    c(B = "Raw B (cell proportion)", B_std = "Standardized B (affinity)"))) +
  labs(x = "Borealness", y = "Number of species") +
  theme_minimal()

idx |> filter(B_std > 0.5)
idx |> filter(B_std < 0.5)

idx |>
  tidyr::pivot_longer(c(B, B_std), names_to = "index", values_to = "value") |>
  dplyr::summarise(
    n_boreal  = sum(value > 0.5),
    n_tundra  = sum(value < 0.5),
    pct_boreal = round(100 * mean(value > 0.5), 1),
    pct_tundra = round(100 * mean(value < 0.5), 1),
    .by = index
  )

ggplot(idx, aes(B, B_std)) +
  geom_abline(linetype = "dashed") +
  geom_point() +
  ggrepel::geom_text_repel(data = \(d) dplyr::filter(d, abs(B - B_std) > 0.1),
                           aes(label = species), size = 2.5) +
  theme_minimal()

## ============================================================
## GREENLAND MAP  (separate relaxed download; index `idx` stays
## the strict circumpolar version — do not mix the two)
## ============================================================

# --- GL download (heavy, run once; no uncertainty filter) ---
# d_gl <- occ_download(
#   pred_in("taxonKey", keys),
#   pred("country", "GL"),
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   pred_in("basisOfRecord", c("PRESERVED_SPECIMEN", "HUMAN_OBSERVATION")),
#   format = "SIMPLE_CSV"
# )
# occ_download_wait(d_gl)
# occ_gl <- occ_download_get(d_gl) |> occ_download_import()
# saveRDS(occ_gl, "data/occ_raw_gl.rds")
# writeLines(as.character(d_gl), "data/download_key_gl.txt")
occ_gl <- readRDS("data/occ_raw_gl.rds")

# --- clean ---
occ_gl_clean <- occ_gl |>
  as.data.frame() |>
  cc_val(lon = "decimalLongitude", lat = "decimalLatitude") |>
  cc_equ(lon = "decimalLongitude", lat = "decimalLatitude") |>
  cc_zero(lon = "decimalLongitude", lat = "decimalLatitude") |>
  cc_cen(lon = "decimalLongitude", lat = "decimalLatitude", buffer = 1000) |>
  cc_inst(lon = "decimalLongitude", lat = "decimalLatitude")

# --- assign cells, join species-level B, summarise per cell ---
pts_gl <- occ_gl_clean |>
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) |>
  st_transform("EPSG:6931") |> vect()

gl_cells <- occ_gl_clean |>
  dplyr::mutate(cell  = cells(biome_r, pts_gl)[, "cell"],
                biome = terra::extract(biome_r, pts_gl)[[2]]) |>
  dplyr::filter(!is.na(biome)) |>
  dplyr::distinct(species, cell) |>
  dplyr::inner_join(idx, by = "species") |>
  dplyr::summarise(mean_B = mean(B_std), n_sp = dplyr::n(), .by = cell)
nrow(gl_cells)   # expect 141

# --- raster fill + plot (one chunk so it can't go stale) ---

gl <- rnaturalearth::ne_countries(country = "Greenland", scale = 50) |>
  st_transform("EPSG:6931")
crs_gl <- "+proj=laea +lat_0=90 +lon_0=-40 +datum=WGS84 +units=m"

b_map <- rast(biome_r); values(b_map) <- NA
b_map[gl_cells$cell] <- gl_cells$mean_B
b_map_gl <- crop(b_map, vect(gl))     # crop only; no mask
bb <- st_bbox(gl_rot)   # bbox in the rotated crs — must come from gl_rot, not gl

ggplot() +
  geom_sf(data = gl_rot, fill = "grey90", color = NA) +
  geom_spatraster(data = b_map_rot) +
  scale_fill_viridis_c(name = "Mean B_std", na.value = NA) +
  coord_sf(crs = crs_gl,
           xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"], bb["ymax"]),
           expand = FALSE) +
  theme_minimal()

