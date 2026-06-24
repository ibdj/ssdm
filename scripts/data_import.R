#### packages ##################################################################

library(tidyverse)
library(vegan)
library(janitor)
library(terra)
library(sf)
library(whitebox) # to calculate twi
whitebox::install_whitebox()  # installs the WhiteboxTools binary
library(readxl)
library(gstat)
#devtools::install_github('cjcarlson/embarcadero') #sdm also for interpolation
library(embarcadero)
library(pROC)
library(dbarts)
library(raster)
library(caret)
library(GGally) # to make correlation plots for the soil moisture
library(spatialEco)

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

# Single processing function for rasters
process_rast <- function(r, ref = ref_rast) {
  r |>
    project("EPSG:32622") |> #common projection
    crop(aoi) |>             #cropping
    resample(ref)            #resampling to the same reference layer
}

#### loading tms data ##########################################################

raw_tms_mp <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/tms_pivot.rds")

tms_mp <- raw_tms_mp |> 
  clean_names() |> 
  filter(level == "t1_below6cm") |>          # soil temperature only
  filter(month(date) %in% c(6, 7, 8)) |>  # June-August         # June-August
  group_by(plot) |> 
  reframe(
    temp_mean_tms = mean(temp, na.rm = TRUE),
    mois_raw_tms = mean(moisture_data_raw, na.rm = TRUE)
  ) |> 
  mutate(
    plot_name = toupper(plot)
  )

summary(tms_mp)

raw_BioBasis_Nuuk_PhenologyPlots_Microclimate_2025 <- read_excel("~/Library/CloudStorage/OneDrive-GrønlandsNaturinstitut/General - BioBasis/03_GEM_Database/Datafiler excel/BioBasis_Nuuk_PhenologyPlots_Microclimate_2025.xlsx")
raw_BioBasis_Nuuk_CFlux_Microclimate_2025 <- read_excel("~/Library/CloudStorage/OneDrive-GrønlandsNaturinstitut/General - BioBasis/03_GEM_Database/Datafiler excel/BioBasis_Nuuk_CFlux_Microclimate_2025.xlsx")

raw_tms_biobasis <- raw_BioBasis_Nuuk_CFlux_Microclimate_2025 |> 
  bind_rows(raw_BioBasis_Nuuk_PhenologyPlots_Microclimate_2025) |> 
  filter(month(Date) %in% c(6, 7, 8)) |>  # June-August
  group_by(Plot, Latitude, Longitude) |>
  reframe(
    temp_mean_tms = mean(Temp_6cmbel, na.rm = TRUE),
    mois_raw_tms = mean(Raw_soil_moisture, na.rm = TRUE)
  )

summary(raw_tms_biobasis)

raw_samples_qgis <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/samples_qgis.csv") |> 
  dplyr::select(plot, X, Y, elevation, ndvi, ndwi) |> 
  mutate(plot_name = plot) |> 
  left_join(tms_mp, by = "plot_name") |> 
  rename(plot = plot.x, plot_tms = plot.y)

names(raw_samples_qgis)

raw_qgis_samples <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/samples.csv", col_types = cols(Date = col_datetime(format = "%m/%d/%Y %H.%M"))) |> 
  clean_names() |> 
  mutate(rowid = row_number(), plot_name = toupper(plot_name))

raw_df_cover <- raw_qgis_samples |> 
 mutate(across(ends_with("_bb"), bb_to_cover)) |> 
  mutate(total_cover = rowSums(across(ends_with("_bb")), na.rm = TRUE))

summary(raw_df_cover)
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

#### abiotic df (field measurements) ################################################################

mp_abiotic <- raw_df_cover |>
  left_join(species_matrix |> dplyr::select(plot_name), by = "plot_name") |>
  mutate(
    richness = rowSums(sp_cols > 0),
    shannon = vegan::diversity(sp_cols, index = "shannon")
  )

summary(mp_abiotic)

mp_abiotic <- mp_abiotic |> 
  dplyr::select(plot_name, veg_height_ave, bare_ground_bb, x, y, total_cover, richness, shannon, soil_moi_ave, soil_tem_ave) |> 
  rename(mois_mean_mea = soil_moi_ave,
         temp_mean_mea = soil_tem_ave) |> 
  left_join(tms_mp, by = "plot_name") |> 
  dplyr::select(-plot)

summary(mp_abiotic)

# Find the plot with NA temp
na_plot <- mp_abiotic |> 
  filter(is.na(temp_mean_mea))

# Convert to sf and extract from raster
na_plot_sf <- na_plot |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(32622)

na_plot_sf

#### combining all tms ##########################################################

# Extract TMS logger plots with coordinates from abiotic_plot
tms_mp_own <-  mp_abiotic |>
  filter(!is.na(temp_mean_tms)) |>
  dplyr::select(plot_name, x, y, temp_mean_tms, mois_raw_tms) |>
  rename(Longitude = x, Latitude = y)

# Bind and normalise rmi together (RELATIVE MOISTURE INDEX)
tms_combined <- bind_rows(tms_mp_own, raw_tms_biobasis |> rename(plot_name = Plot)) |>
  mutate(
    rmi_tms = (mois_raw_tms - min(mois_raw_tms)) / 
      (max(mois_raw_tms) - min(mois_raw_tms)) * 100
  )

summary(tms_combined)

nrow(tms_combined)

tms_combined_sf <- tms_combined |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) |>
  st_transform(32622)

#### species frequency #########################################################

species_frequency <- species_matrix |>
  dplyr::select(-plot_name) |>
  summarise(across(everything(), ~ sum(. > 0))) |>
  pivot_longer(everything(), 
               names_to = "species", 
               values_to = "n_plots") |>
  arrange(n_plots)

print(species_frequency, n = Inf)

#### plots and aoi #############################################################

plots_sf <- mp_abiotic |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(32622)

aoi <- plots_sf |>
  st_bbox() |>
  st_as_sfc() |>
  st_buffer(50) |>
  vect()  # convert to terra format for cropping

#### aoi export to python/gee ####################

# aoi_sf <- plots_sf |>
#   st_bbox() |>
#   st_as_sfc() |>
#   st_buffer(50) |>
#   st_transform(4326)  # GEE needs WGS84
# 
# st_write(aoi_sf, "data/aoi.shp", delete_dsn = TRUE)
# #This creates 4 files (.shp, .shx, .dbf, .prj) — you need to upload all four to GEE as a zip:
#   # Zip all shapefile components for GEE upload
# zip("data/aoi.zip", 
#     files = c("data/aoi.shp", "data/aoi.shx", 
#               "data/aoi.dbf", "data/aoi.prj"))

#### raster import #############################################################

dem_rast        <- rast("data/elevation_arcticdem-30_32622.tif")
ndvi_rast       <- rast("data/ndvi_export_2025.tif")
ndwi_rast       <- rast("data/ndwi.tif")
snowfree_rast   <- rast("data/snow_free_days.tif")

slope_rast      <- terrain(dem_rast, v = "slope", unit = "degrees")
aspect_rast     <- terrain(dem_rast, v = "aspect", unit = "degrees")
aspect_cos_rast <- cos(aspect_rast * pi / 180)
aspect_sin_rast <- sin(aspect_rast * pi / 180)

summary(ndwi_rast)

#### raster solar radiation ####################################################

#### raster standardising 1 ####################################################

# Define reference raster - everything gets matched to this
ref_rast <- rast("data/ndvi_export_2025.tif") |>
  project("EPSG:32622") |>
  crop(aoi)

rast_dem_proc        <- dem_rast |> process_rast()
rast_ndvi_proc       <- ndvi_rast |> process_rast()
rast_ndwi_proc       <- ndwi_rast |> process_rast()
rast_snowfree_proc   <- snowfree_rast |> process_rast()
rast_slope_proc      <- slope_rast |> process_rast()
rast_aspect_proc     <- aspect_rast |> process_rast()
rast_aspect_cos_proc <- aspect_cos_rast |> process_rast()
rast_aspect_sin_proc <- aspect_sin_rast |> process_rast()

sapply(list(rast_dem_proc, 
            rast_ndvi_proc, 
            rast_ndwi_proc, 
            rast_snowfree_proc,
            rast_slope_proc,
            rast_aspect_proc,
            rast_aspect_cos_proc,
            rast_aspect_sin_proc
            ),
function(r) crs(r, describe = TRUE)$code)

#### raster twi (calculating) ##################################################

#commented out because I will use ndwi instead

# wbt_fill_depressions("data/dem_crop.tif", "data/dem_filled.tif")
# wbt_d8_flow_accumulation("data/dem_filled.tif", "data/sca.tif")
# wbt_slope("data/dem_filled.tif", "data/slope_wb.tif", units = "degrees")
# wbt_wetness_index(
#   sca = "data/sca.tif",
#   slope = "data/slope_wb.tif",
#   output = "data/twi_calculated.tif"
# )
# 
# twi_rast <- rast("data/twi_calculated.tif")

#### sampling all imported rasters #############################################

tms_combined <- tms_combined |>
  mutate(
    elevation  = terra::extract(rast_dem_proc,        tms_combined_sf)[, 2],
    ndvi       = terra::extract(rast_ndvi_proc,       tms_combined_sf)[, 2],
    ndwi       = terra::extract(rast_ndwi_proc,       tms_combined_sf)[, 2],
    snowfree   = terra::extract(rast_snowfree_proc,   tms_combined_sf)[, 2],
    slope      = terra::extract(rast_slope_proc,      tms_combined_sf)[, 2],
    aspect_raw = terra::extract(rast_aspect_proc,     tms_combined_sf)[, 2],
    aspect_cos = terra::extract(rast_aspect_cos_proc, tms_combined_sf)[, 2],
    aspect_sin = terra::extract(rast_aspect_sin_proc, tms_combined_sf)[, 2],
  )

summary(tms_combined)

#### raster temp (interpolation) ###############################################

# Step 1: fit linear model on combined logger data
temp_lm <- lm(temp_mean_tms ~ ndvi + aspect_cos + aspect_sin + elevation, 
              data = tms_combined)

summary(temp_lm)

# Standard linear model diagnostic plots
par(mfrow = c(2, 2))
plot(temp_lm)
par(mfrow = c(1, 1))

# Add residuals to combined logger data
tms_combined_sf <- tms_combined_sf |>
  mutate(temp_resid = residuals(temp_lm))

# Compute variogram of residuals
vgm_temp <- variogram(temp_resid ~ 1, data = tms_combined_sf)
plot(vgm_temp)

# Now stack
pred_stack <- c(rast_ndvi_proc, rast_aspect_cos_proc, rast_aspect_sin_proc, rast_dem_proc)
names(pred_stack) <- c("ndvi", "aspect_cos", "aspect_sin", "elevation")

# Project
temp_rast <- predict(pred_stack, temp_lm)
plot(temp_rast)

#### raster moisture (just checking the bad correlation) #######################

my_scatter <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, 
                colour = "red", linewidth = 0.8)
}

tms_combined |> 
  dplyr::select(
    `RMI (logger)` = rmi_tms,
    `Soil moisture (field)` = mois_mean_mea,
    `NDWI` = ndwi
  ) |>
  ggpairs(
    upper = list(continuous = wrap("cor", size = 4)),
    lower = list(continuous = my_scatter),
    diag  = list(continuous = wrap("densityDiag"))
  ) +
  labs(title = "Correlation between moisture measures") +
  theme_minimal()

#### raster processing 2 #######################################################

rast_twi_proc        <- twi_rast |> process_rast()
rast_temp_proc       <- temp_rast |> process_rast()

#### checking the nas ##########################################################

abiotic_plot |> 
  filter(is.na(elevation) | is.na(slope) | is.na(aspect_raw)) |> 
  dplyr::select(plot_name, x, y, elevation, slope, aspect_raw, twi, ndvi)

abiotic_plot |> 
  filter(is.na())

summary(abiotic_plot)


##### interpolation ############################################################

tms_combined_sf <- tms_combined |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) |>
  st_transform(32622)

tms_combined <- tms_combined |>
  mutate(
    elevation = terra::extract(dem_rast, tms_combined_sf)[, 2],
    slope = terra::extract(slope_rast, tms_combined_sf)[, 2],
    aspect_raw = terra::extract(aspect_rast, tms_combined_sf)[, 2],
    aspect_sin = sin(aspect_raw * pi / 180),
    aspect_cos = cos(aspect_raw * pi / 180),
    twi = terra::extract(twi_rast, tms_combined_sf)[, 2],
    ndvi = terra::extract(ndvi_rast, tms_combined_sf)[, 2]
  )

# Recheck correlations
tms_combined |>
  dplyr::select(temp_mean_tms, vwc_tms, elevation, slope, twi, ndvi, aspect_sin, aspect_cos) |>
  cor(use = "complete.obs") |>
  round(2)


#testing for correlation
abiotic_plot |>
  dplyr::select(soil_tem_ave, soil_moi_ave, elevation, slope, twi, ndvi, aspect_sin, aspect_cos) |>
  cor(use = "complete.obs") |>
  round(2)

abiotic_plot |>
  dplyr::select(temp_mean_tms, vwc_tms, elevation, slope, twi, ndvi, aspect_sin, aspect_cos) |>
  cor(use = "complete.obs") |>
  round(2)

#This will tell if the logger-based values show cleaner relationships with topography than the point measurements.


#### Extract interpolated value ################################################
imputed_temp <- terra::extract(temp_rast, na_plot_sf)[, 2]
imputed_temp

# Fill NA in abiotic_plot
abiotic_plot <- abiotic_plot |>
  mutate(temp_predicted = ifelse(is.na(soil_tem_ave), 
                                 imputed_temp, 
                                 soil_tem_ave))

# Verify no more NAs
sum(is.na(abiotic_plot$temp_predicted))

#### raster sampling ###########################################################

abiotic_plot <- abiotic_plot |>
  mutate(
    elevation  = terra::extract(rast_dem_proc,        plots_sf)[, 2],
    ndvi       = terra::extract(rast_ndvi_proc,       plots_sf)[, 2],
    ndwi       = terra::extract(rast_ndwi_proc,       plots_sf)[, 2],
    snowfree   = terra::extract(rast_snowfree_proc,   plots_sf)[, 2],
    slope      = terra::extract(rast_slope_proc,      plots_sf)[, 2],
    aspect_raw = terra::extract(rast_aspect_proc,     plots_sf)[, 2],
    aspect_cos = terra::extract(rast_aspect_cos_proc, plots_sf)[, 2],
    aspect_sin = terra::extract(rast_aspect_sin_proc, plots_sf)[, 2],
    twi        = terra::extract(rast_twi_proc,        plots_sf)[, 2],
    temp       = terra::extract(rast_temp_proc,       plots_sf)[, 2]
  )
#### writing all data files ####################################################
saveRDS(abiotic_plot, "data/abiotic_plot.rds")
saveRDS(species_matrix, "data/species_matrix.rds")
saveRDS(species_long, "data/species_long.rds")
write_rds(species_frequency, "data/species_frequency.rds")

writeRaster(rast_dem_proc, "data/rast_dem_proc.tif", overwrite = TRUE)
writeRaster(rast_ndvi_proc, "data/rast_ndvi_proc.tif", overwrite = TRUE)
writeRaster(rast_ndwi_proc, "data/rast_ndwi_proc.tif", overwrite = TRUE)
writeRaster(rast_snowfree_proc, "data/rast_snowfree_proc.tif", overwrite = TRUE)
writeRaster(rast_slope_proc, "data/rast_slope_proc.tif", overwrite = TRUE)
writeRaster(rast_aspect_proc, "data/rast_aspect_proc.tif", overwrite = TRUE)
writeRaster(rast_aspect_cos_proc, "data/rast_aspect_cos_proc.tif", overwrite = TRUE)
writeRaster(rast_aspect_sin_proc, "data/rast_aspect_sin_proc.tif", overwrite = TRUE)
writeRaster(rast_twi_proc, "data/rast_twi_proc.tif", overwrite = TRUE)
