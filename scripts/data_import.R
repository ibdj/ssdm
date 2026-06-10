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

# Single processing function
process_rast <- function(r, ref = ref_rast) {
  r |>
    project("EPSG:32622") |>
    crop(aoi) |>
    resample(ref)
}

#### loading tms data ##########################################################

tms_mp <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/tms_pivot.rds") |> 
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

BioBasis_Nuuk_PhenologyPlots_Microclimate_2025 <- read_excel("~/Library/CloudStorage/OneDrive-GrønlandsNaturinstitut/General - BioBasis/03_GEM_Database/Datafiler excel/BioBasis_Nuuk_PhenologyPlots_Microclimate_2025.xlsx")
BioBasis_Nuuk_CFlux_Microclimate_2025 <- read_excel("~/Library/CloudStorage/OneDrive-GrønlandsNaturinstitut/General - BioBasis/03_GEM_Database/Datafiler excel/BioBasis_Nuuk_CFlux_Microclimate_2025.xlsx")

tms_biobasis <- BioBasis_Nuuk_CFlux_Microclimate_2025 |> 
  bind_rows(BioBasis_Nuuk_PhenologyPlots_Microclimate_2025) |> 
  filter(month(Date) %in% c(6, 7, 8)) |>  # June-August
  group_by(Plot, Latitude, Longitude) |>
  reframe(
    temp_mean_tms = mean(Temp_6cmbel, na.rm = TRUE),
    mois_raw_tms = mean(Raw_soil_moisture, na.rm = TRUE)
  )

summary(tms_biobasis)

samples_qgis <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/samples_qgis.csv") |> 
  dplyr::select(plot, X, Y, elevation, ndvi, ndwi) |> 
  mutate(plot_name = plot) |> 
  left_join(tms_mp, by = "plot_name") |> 
  rename(plot = plot.x, plot_tms = plot.y)

names(samples_qgis)

df_raw <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/samples.csv", col_types = cols(Date = col_datetime(format = "%m/%d/%Y %H.%M"))) |> 
  clean_names() |> 
  mutate(rowid = row_number(), plot_name = toupper(plot_name))

df_cover <- df_raw |> 
 mutate(across(ends_with("_bb"), bb_to_cover)) |> 
  mutate(total_cover = rowSums(across(ends_with("_bb")), na.rm = TRUE))

summary(df_cover)
#### species matrix ############################################################

# Step 1: pivot just the species names to long
taxon_names <- df_cover |>
  dplyr::select(plot_name, matches("^taxon_[0-9]+$")) |>
  pivot_longer(-plot_name, names_to = "slot", values_to = "species_name")

# Step 2: pivot just the bb values to long
taxon_bb <- df_cover |>
  dplyr::select(plot_name, matches("^taxon_[0-9]+_bb$")) |>
  pivot_longer(-plot_name, names_to = "slot", values_to = "cover") |>
  mutate(slot = str_remove(slot, "_bb$"))

# Step 3: pivot just the height values to long
taxon_height <- df_cover |>
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

#### abiotic df ################################################################

abiotic_plot <- df_cover |>
  left_join(species_matrix |> dplyr::select(plot_name), by = "plot_name") |>
  mutate(
    richness = rowSums(sp_cols > 0),
    shannon = vegan::diversity(sp_cols, index = "shannon")
  )

abiotic_plot <- abiotic_plot |> 
  dplyr::select(plot_name, veg_height_ave, bare_ground_bb, x, y, total_cover, richness, shannon, soil_moi_ave, soil_tem_ave) |> 
  rename(mois_mean_mea = soil_moi_ave,
         temp_mean_mea = soil_tem_ave) |> 
  left_join(tms_mp, by = "plot_name") |> 
  dplyr::select(-plot)

summary(abiotic_plot)

# Find the plot with NA temp
na_plot <- abiotic_plot |> 
  filter(is.na(temp_mean_mea))

# Convert to sf and extract from raster
na_plot_sf <- na_plot |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(32622)

na_plot_sf

#### combining all tms ##########################################################

# Extract TMS logger plots with coordinates from abiotic_plot
tms_own <-  abiotic_plot |>
  filter(!is.na(temp_mean_tms)) |>
  dplyr::select(plot_name, x, y, temp_mean_tms, mois_raw_tms) |>
  rename(Longitude = x, Latitude = y)

# Bind and normalise rmi together (RELATIVE MOISTURE INDEX)
tms_combined <- bind_rows(tms_own, tms_biobasis) |>
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

write_rds(species_frequency, "data/species_frequency.rds")


#### plots and aoi #############################################################

plots_sf <- abiotic_plot |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(32622)

aoi <- plots_sf |>
  st_bbox() |>
  st_as_sfc() |>
  st_buffer(50) |>
  vect()  # convert to terra format for cropping

#### raster imports ############################################################

dem_rast        <- rast("data/elevation_arcticdem-30_32622.tif")
ndvi_rast       <- rast("data/ndvi_export_2025.tif")
ndwi_rast       <- rast("data/ndwi.tif")
snowfree_rast   <- rast("data/snow_free_days.tif")
slope_rast      <- terrain(dem_rast, v = "slope", unit = "degrees")
aspect_rast     <- terrain(dem_rast, v = "aspect", unit = "degrees")
aspect_cos_rast <- cos(aspect_rast * pi / 180)
aspect_sin_rast <- sin(aspect_rast * pi / 180)
  
#### raster twi (calculating) ##################################################

wbt_fill_depressions("data/dem_crop.tif", "data/dem_filled.tif")
wbt_d8_flow_accumulation("data/dem_filled.tif", "data/sca.tif")
wbt_slope("data/dem_filled.tif", "data/slope_wb.tif", units = "degrees")
wbt_wetness_index(
  sca = "data/sca.tif",
  slope = "data/slope_wb.tif",
  output = "data/twi_calculated.tif"
)

twi_rast <- rast("data/twi_calculated.tif")

#### raster temp (interpolation) ###############################################

# Step 1: fit linear model on combined logger data
temp_lm <- lm(temp_mean_tms ~ ndvi + aspect_cos + aspect_sin, 
              data = tms_combined)

summary(temp_lm)

# Add residuals to combined logger data
tms_combined <- tms_combined |>
  mutate(temp_resid = residuals(temp_lm))

# Convert to sf for variogram
tms_combined_sf <- tms_combined |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) |>
  st_transform(32622)

# Compute variogram of residuals
vgm_temp <- variogram(temp_resid ~ 1, data = tms_combined_sf)
plot(vgm_temp)

# Now stack
pred_stack <- c(rast_ndvi_proc, rast_aspect_cos_proc, rast_aspect_sin_proc)
names(pred_stack) <- c("ndvi", "aspect_cos", "aspect_sin")

# Project
temp_rast <- predict(pred_stack, temp_lm)
plot(temp_rast)

#temp_rast_masked <- mask(temp_rast, ndvi_rast < 0.1, maskvalue = TRUE)
#plot(temp_rast_masked)
#writeRaster(temp_rast_masked, "data/temp_predicted_rast.tif", overwrite = TRUE)

#### moisture interpolation ####################################################

moisture_sf <- abiotic_plot |>
  filter(!is.na(soil_moi_ave)) |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(32622)

vgm_moist <- variogram(soil_moi_ave ~ 1, data = moisture_sf)
plot(vgm_moist)
vgm_moist_fit <- fit.variogram(vgm_moist, 
                               model = vgm(psill = 300, 
                                           model = "Sph", 
                                           range = 800, 
                                           nugget = 150))
plot(vgm_moist, vgm_moist_fit)

# Create prediction grid from raster extent
pred_grid <- as.points(ndvi_rast) |>
  st_as_sf() |>
  st_transform(32622)

# Run ordinary kriging
moist_krige <- krige(
  formula = soil_moi_ave ~ 1,
  locations = moisture_sf,
  newdata = pred_grid,
  model = vgm_moist_fit
)

# Check what moist_krige looks like
class(moist_krige)
head(moist_krige)
nrow(moist_krige)
summary(moist_krige$var1.pred)

# Extract coordinates and predictions
moist_df <- st_coordinates(moist_krige) |>
  as.data.frame() |>
  mutate(moisture = moist_krige$var1.pred)

# Convert to raster
moist_rast <- rast(moist_df, type = "xyz", crs = crs(ndvi_rast))

# # Check
print(moist_rast)
summary(moist_rast)
plot(moist_rast)

# Mask water bodies same as temperature
moist_rast_masked <- mask(moist_rast, ndvi_rast < 0.1, maskvalue = TRUE)

plot(moist_rast_masked)

# Simple regression on NDVI first
moist_lm4 <- lm(soil_moi_ave ~ ndvi + twi + elevation, data = abiotic_plot)
summary(moist_lm4)

# Prepare data - use all available predictors
moist_bart <- bart(
  x.train = abiotic_plot |> 
    dplyr::select(ndvi, twi, elevation, slope, aspect_sin, aspect_cos, snowfree) |>
    as.data.frame(),
  y.train = abiotic_plot$soil_moi_ave,
  keeptrees = TRUE
)

summary(moist_bart)


#### processing all rasters ####################################################

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
rast_twi_proc        <- twi_rast |> process_rast()
rast_temp_proc       <- temp_rast |> process_rast()

sapply(list(rast_dem_proc, 
            rast_ndvi_proc, 
            rast_ndwi_proc, 
            rast_snowfree_proc,
            rast_slope_proc,
            rast_aspect_proc,
            rast_aspect_cos_proc,
            rast_aspect_sin_proc,
            rast_twi_proc,
            rast_temp_proc
            ),
       function(r) crs(r, describe = TRUE)$code)

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


#### writing all data files ####################################################
saveRDS(abiotic_plot, "data/abiotic_plot.rds")
saveRDS(species_matrix, "data/species_matrix.rds")
saveRDS(species_long, "data/species_long.rds")

writeRaster(rast_dem_proc, "data/rast_dem_proc.tif", overwrite = TRUE)
writeRaster(rast_ndvi_proc, "data/rast_ndvi_proc.tif", overwrite = TRUE)
writeRaster(rast_ndwi_proc, "data/rast_ndwi_proc.tif", overwrite = TRUE)
writeRaster(rast_snowfree_proc, "data/rast_snowfree_proc.tif", overwrite = TRUE)
writeRaster(rast_slope_proc, "data/rast_slope_proc.tif", overwrite = TRUE)
writeRaster(rast_aspect_proc, "data/rast_aspect_proc.tif", overwrite = TRUE)
writeRaster(rast_aspect_cos_proc, "data/rast_aspect_cos_proc.tif", overwrite = TRUE)
writeRaster(rast_aspect_sin_proc, "data/rast_aspect_sin_proc.tif", overwrite = TRUE)
writeRaster(rast_twi_proc, "data/rast_twi_proc.tif", overwrite = TRUE)
