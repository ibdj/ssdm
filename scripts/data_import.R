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
library(spatialEco) # for radiation / heat load index
library(car)
library(ranger) # to do regression trees on the soil moisture from tms
library(gstat)

#### functions #################################################################


#cover function needed to convert the bryophyte and lichen cover
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


##### interpolation the soil moisture ##########################################

vwc_agg <- tms_long_filt |>
  dplyr::filter(sensor_name == "VWC_universal") |>
  dplyr::filter(lubridate::month(datetime) %in% 6:8) |>
  dplyr::filter(!is.na(value)) |>
  dplyr::summarise(
    vwc = mean(value, na.rm = TRUE),
    n_obs = dplyr::n(),
    .by = c(serial, plot)
  )

summary(vwc_agg)
vwc_agg 

################################################################################

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
#### abiotic df (field measurements) ################################################################

names(survey_0_ren)
plot_meta <- survey_0_ren |>
  dplyr::select(1:24) |>
  dplyr::rename(bryophyte_bb = bryophyte_bb_49) |>
  dplyr::mutate(
    cover_bryophyte = bb_to_cover(bryophyte_bb),
    cover_lichen    = bb_to_cover(lichen_bb),
    cove_bareground = bb_to_cover(bare_ground_bb)
  )

plot_meta <- plot_meta |>
  dplyr::mutate(dplyr::across(
    dplyr::matches("^soil_(moisture|temp)_[news]$"),
    ~ dplyr::if_else(.x < 0, NA_real_, .x)
  ))

plot_meta <- plot_meta |>
  dplyr::mutate(veg_height_mean = rowMeans(
    dplyr::pick(dplyr::starts_with("vegetation_height_")), na.rm = TRUE
  ))

# calculating the veg height from the species data
veg_species_height <- species_long |>
  dplyr::group_by(plot_name) |>
  dplyr::summarise(cover_allplant = round(sum(cover, na.rm = TRUE),2), veg_height_species = round(mean(height, na.rm = TRUE),2), ) |>
  dplyr::right_join(dplyr::distinct(survey_0, plot_name), by = "plot_name") |>
  dplyr::mutate(dplyr::across(-plot_name, ~ tidyr::replace_na(.x, 0)))

cover_functional_groups <- species_long |>
  dplyr::group_by(plot_name, func_type) |>
  dplyr::summarise(cover = sum(cover, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(
    names_from = func_type,
    values_from = cover,
    names_prefix = "cover_",
    values_fill = 0
  ) 
  
# joining species mean height and the plot veg height
plot_meta <- plot_meta |> 
  left_join(veg_species_height, by = "plot_name") |> 
  left_join(cover_functional_groups, by = "plot_name") |> 
  dplyr::select(plot_name,x,y,rowid,cover_bryophyte,cover_lichen,cove_bareground,veg_height_mean,veg_height_species, cover_allplant, cover_forb, cover_graminoid, cover_shrub)

# Convert to sf and extract from raster
plot_meta <- plot_meta |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(32622)

summary(plot_meta)

#basic plotting
ggplot(plot_meta, aes(x = cover_allplant, y = veg_height_species)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Vegetation height, plot corners (cm)",
    y = "Vegetation height, species mean (cm)"
  ) +
  theme_minimal()


#### species frequency #########################################################

species_frequency <- species_matrix |>
  dplyr::select(-plot_name) |>
  summarise(across(everything(), ~ sum(. > 0))) |>
  pivot_longer(everything(), 
               names_to = "species", 
               values_to = "n_plots") |>
  arrange(n_plots)

print(species_frequency, n = Inf)


#plot(aoi_buffered, border = "blue")
#plot(aoi_masked, add = TRUE, border = "red")
#plot(plots_sf, add = TRUE)
#### aoi export to python/gee ##################################################


#### raster standardising 1 ####################################################

# Define reference raster - everything gets matched to this
ref_rast <- rast("data/ndvi_export_2025.tif") |>
  mask(aoi_masked) |> 
  project("EPSG:32622")
  
plot(trim(ref_rast))
plot(aoi_masked, add = TRUE)

rast_dem_proc        <- rast_dem |> process_rast()
rast_ndvi_proc       <- rast_ndvi |> process_rast()
rast_ndwi_proc       <- rast_ndwi |> process_rast()
rast_snowfree_proc   <- rast_snowfree |> process_rast()
rast_slope_proc      <- rast_slope |> process_rast()
#rast_aspect_proc     <- rast_aspect |> process_rast()
#rast_aspect_cos_proc <- rast_aspect_cos |> process_rast()
#rast_aspect_sin_proc <- rast_aspect_sin |> process_rast()

#rast_tpi_proc        <- tpi |> process_rast()
rast_hli_proc        <- hli |> process_rast()
rast_temp_proc       <- temp_rast |> process_rast()

sapply(list(rast_dem_proc, 
            rast_ndvi_proc, 
            rast_ndwi_proc, 
            rast_snowfree_proc,
            rast_slope_proc,
            rast_aspect_proc,
            rast_aspect_cos_proc,
            rast_aspect_sin_proc,
            rast_tpi_proc,
            rast_hli_proc,
            rast_temp_proc
            ),
function(r) crs(r, describe = TRUE)$code)

plot(trim(rast_snowfree_proc))
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

#### sampling based on the new all tms file ####################################

tms_sf <- tms_sf |>
  mutate(
    elevation = terra::extract(rast_dem_proc,      tms_sf)[, 2],
    ndvi      = terra::extract(rast_ndvi_proc,     tms_sf)[, 2],
    ndwi      = terra::extract(rast_ndwi_proc,     tms_sf)[, 2],
    snowfree  = terra::extract(rast_snowfree_proc, tms_sf)[, 2],
    slope     = terra::extract(rast_slope_proc,    tms_sf)[, 2],
#    tpi       = terra::extract(rast_tpi_proc,      tms_sf)[, 2],
    hli       = terra::extract(rast_hli_proc,      tms_sf)[, 2]
  )

vwc_model <- tms_sf |>
  mutate(serial = as.character(serial)) |>
  dplyr::left_join(vwc_agg |> mutate(serial = as.character(serial)),
                   by = "serial") |>
  sf::st_drop_geometry()

colSums(is.na(vwc_model[, c("vwc","elevation","ndvi","ndwi","snowfree","slope","hli")]))

#### explaingning the vwc by variables ####

cand <- c("elevation", "ndvi", "ndwi", "snowfree", "slope", "hli")

# 1. how does each predictor relate to VWC?
round(cor(vwc_model[, c("vwc", cand)]), 2)

ctrl <- trainControl(method = "LOOCV", allowParallel = FALSE)

full_vwc <- train(vwc ~ elevation + ndvi + ndwi + snowfree + slope + hli,
                  data = vwc_model, method = "lm", trControl = ctrl)
full_vwc$results   # RMSE + R²

# baseline: just predicting the mean
sd(vwc_model$vwc)  # if full-model RMSE ≈ this, the model adds nothing

# regression tree to see non lniear relations ship

library(caret)
# install.packages("ranger")  # if not already installed

ctrl <- trainControl(method = "LOOCV", allowParallel = FALSE)

set.seed(42)
rf_vwc <- train(
  vwc ~ elevation + ndvi + ndwi + snowfree + slope + hli,
  data = vwc_model,
  method = "ranger",
  trControl = ctrl,
  tuneLength = 3,            # tries a few mtry/splitrule settings
  num.trees = 500
)

rf_vwc$results
rf_vwc$bestTune              # which settings won

#### sampling all imported rasters #############################################

tms_combined <- tms_combined |>
  mutate(
    elevation  = terra::extract(rast_dem_proc,        tms_combined_sf)[, 2],
    ndvi       = terra::extract(rast_ndvi_proc,       tms_combined_sf)[, 2],
    ndwi       = terra::extract(rast_ndwi_proc,       tms_combined_sf)[, 2],
    snowfree   = terra::extract(rast_snowfree_proc,   tms_combined_sf)[, 2],
    slope      = terra::extract(rast_slope_proc,      tms_combined_sf)[, 2],
#    aspect_raw = terra::extract(rast_aspect_proc,     tms_combined_sf)[, 2],
#    aspect_cos = terra::extract(rast_aspect_cos_proc, tms_combined_sf)[, 2],
#    aspect_sin = terra::extract(rast_aspect_sin_proc, tms_combined_sf)[, 2],
#    tpi        = terra::extract(rast_tpi_proc,        tms_combined_sf)[, 2],
    hli        = terra::extract(rast_hli_proc,        tms_combined_sf)[, 2]
  )

summary(tms_combined)

#### raster temp (interpolation) ###############################################

# --- Temperature interpolation: predictor selection ---
# Goal: pick a small, physically sensible set of spatial predictors to
# interpolate mean growing-season soil temperature across the AOI.

# --- Aggregate soil temperature per logger (growing season mean) ---
temp_agg <- tms_long_filt |>
  dplyr::filter(sensor_name == "TMS_T1") |>              # soil temp, 8 cm depth
  dplyr::filter(lubridate::month(datetime) %in% 6:8) |>  # June–August
  dplyr::summarise(
    temp  = mean(value, na.rm = TRUE),
    n_obs = dplyr::n(),
    .by = c(serial, plot)
  )

summary(temp_agg)

# --- Join temperature response to the spatial-predictor object ---
tms_sf <- tms_sf |>
  dplyr::mutate(serial = as.character(serial)) |>
  dplyr::left_join(
    temp_agg |> dplyr::select(serial, temp),
    by = "serial"
  )

sum(is.na(tms_sf$temp))   # should be 0

# check it landed and is complete
sum(is.na(tms_sf$temp))   # should be 0
names(tms_sf)

#checking for multicolenearity
cand <- c("elevation", "hli", "ndvi", "ndwi", "snowfree", "slope")

# Correlation matrix among candidates: flags redundant predictors
# (e.g. ndvi/ndwi). Informs interpretation, not inclusion by itself.
round(cor(sf::st_drop_geometry(tms_sf)[, cand]), 2)

# VIF on the full candidate model: quantifies multicollinearity.
# Rule of thumb: >5 worth a look, >10 severe. Low VIF = safe to include.
full <- lm(temp ~ elevation + hli + ndvi + ndwi + snowfree + slope,
           data = tms_sf)
vif(full)

# --- Best-subset selection by cross-validated prediction error ---
# Fit every non-empty subset of candidates; rank by LOOCV-RMSE.
# Out-of-sample error (not in-sample R2) is the selection criterion.

ctrl <- trainControl(method = "LOOCV", allowParallel = FALSE)

preds <- c("elevation", "hli", "ndvi", "ndwi", "snowfree", "slope")

# all non-empty subsets of the candidate predictors
combos <- unlist(
  lapply(seq_along(preds), \(k) combn(preds, k, simplify = FALSE)),
  recursive = FALSE
)

# LOOCV-RMSE for each candidate model
results <- sapply(combos, function(vars) {
  f <- reformulate(vars, response = "temp")
  train(f, data = tms_sf, method = "lm", trControl = ctrl)$results$RMSE
})

out <- data.frame(
  vars = sapply(combos, paste, collapse = " + "),
  n    = sapply(combos, length),
  RMSE = results
)

out[order(out$RMSE), ] |> head(10)   # 10 best models by CV-RMSE
# making a plot to visualise the lowest RMSE

# --- Fit the selected temperature model ---
# elevation + hli + ndvi + ndwi: lowest LOOCV-RMSE, all physically motivated,
# comfortable for n = 69.
temp_lm <- lm(temp ~ elevation + hli + ndvi + ndwi, data = tms_sf)
summary(temp_lm)

# Diagnostic plots: check linearity, normality, homoscedasticity, influence.
par(mfrow = c(2, 2))
plot(temp_lm)
par(mfrow = c(1, 1))

#### --- Residual variogram: is there spatial structure left for kriging? --- ####

# attach residuals to the spatial object
tms_sf$temp_resid <- residuals(temp_lm)

# empirical variogram of residuals
vgm_temp <- variogram(temp_resid ~ 1, data = tms_sf)
plot(vgm_temp)

# --- Predict temperature surface across the AOI ---
# Residual variogram showed no spatial structure (pure nugget), so plain
# regression prediction is appropriate — no kriging needed.

# stack predictors; names MUST match the model's variables
pred_stack <- c(rast_dem_proc, rast_hli_proc, rast_ndvi_proc, rast_ndwi_proc)
names(pred_stack) <- c("elevation", "hli", "ndvi", "ndwi")

# verify names/ranges before predicting (guards against layer-order mixups)
pred_stack

temp_rast <- terra::predict(pred_stack, temp_lm)
plot(trim(temp_rast))

# --- Finalise and export the temperature layer ---
# Align to reference grid + mask to AOI, matching the other predictors.
rast_temp_proc <- temp_rast |> process_rast()
plot(trim(rast_temp_proc))

# write to disk as GeoTIFF
terra::writeRaster(
  rast_temp_proc,
  filename = "data/rast_temp_proc.tif",   # adjust path as needed
  overwrite = TRUE
)

#### raster interpolation of GROWING SEASON LENGTH #############################

gsl_doy <- tms_long_filt |>
  dplyr::filter(sensor_name == "TMS_T1") |>
  dplyr::mutate(doy = yday(datetime)) |>
  # mean soil temp per logger per day-of-year (pooled across years)
  dplyr::summarise(doy_mean = mean(value, na.rm = TRUE),
                   .by = c(serial, plot, doy)) |>
  # GSL = number of DOY with pooled mean > 0
  dplyr::summarise(gsl = sum(doy_mean > 0, na.rm = TRUE),
                   n_doy = dplyr::n(),      # how many DOY had any data (coverage check)
                   .by = c(serial, plot))

summary(gsl_doy)
gsl_doy

# --- Join GSL response to the spatial-predictor object ---
tms_sf <- tms_sf |>
  dplyr::left_join(gsl_doy |> dplyr::select(serial, gsl), by = "serial")

sum(is.na(tms_sf$gsl))   # should be 0

cand <- c("elevation", "hli", "ndvi", "ndwi", "snowfree", "slope")

round(cor(sf::st_drop_geometry(tms_sf)[, c("gsl", cand)]), 2)

full_gsl <- lm(gsl ~ elevation + hli + ndvi + ndwi + snowfree + slope, data = tms_sf)
vif(full_gsl)

# finding the model for gsl interpolation

ctrl <- trainControl(method = "LOOCV", allowParallel = FALSE)

preds <- c("elevation", "hli", "ndvi", "ndwi", "snowfree", "slope")

combos <- unlist(
  lapply(seq_along(preds), \(k) combn(preds, k, simplify = FALSE)),
  recursive = FALSE
)

results <- sapply(combos, function(vars) {
  f <- reformulate(vars, response = "gsl")
  train(f, data = tms_sf, method = "lm", trControl = ctrl)$results$RMSE
})

out_gsl <- data.frame(
  vars = sapply(combos, paste, collapse = " + "),
  n    = sapply(combos, length),
  RMSE = results
)

out_gsl[order(out_gsl$RMSE), ] |> head(10)

sd(tms_sf$gsl)

# R² of the best model
best_gsl <- train(gsl ~ hli + ndvi, data = tms_sf,
                  method = "lm", trControl = ctrl)
best_gsl$results   # look at Rsquared

set.seed(42)
rf_gsl <- train(
  gsl ~ elevation + hli + ndvi + ndwi + snowfree + slope,
  data = tms_sf, method = "ranger",
  trControl = ctrl, tuneLength = 3, num.trees = 500
)
rf_gsl$results

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

#### checking the nas ##########################################################

mp_abiotic |> 
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

#### raster sampling ###########################################################

mp_abiotic <- mp_abiotic |>
  mutate(
    elevation  = terra::extract(rast_dem_proc,        plots_sf)[, 2],
    ndvi       = terra::extract(rast_ndvi_proc,       plots_sf)[, 2],
    ndwi       = terra::extract(rast_ndwi_proc,       plots_sf)[, 2],
    snowfree   = terra::extract(rast_snowfree_proc,   plots_sf)[, 2],
    slope      = terra::extract(rast_slope_proc,      plots_sf)[, 2],
    aspect_raw = terra::extract(rast_aspect_proc,     plots_sf)[, 2],
    aspect_cos = terra::extract(rast_aspect_cos_proc, plots_sf)[, 2],
    aspect_sin = terra::extract(rast_aspect_sin_proc, plots_sf)[, 2],
    #twi        = terra::extract(rast_twi_proc,        plots_sf)[, 2],
    temp       = terra::extract(rast_temp_proc,       plots_sf)[, 2],
    tpi        = terra::extract(rast_tpi_proc,        plots_sf)[, 2],
    hli        = terra::extract(rast_hli_proc,       plots_sf)[, 2]
  )
#### writing all data files ####################################################
saveRDS(mp_abiotic, "data/abiotic_plot.rds")
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
writeRaster(rast_hli_proc, "data/rast_hli_proc.tif", overwrite = TRUE)
writeRaster(rast_tpi_proc, "data/rast_tpi_proc.tif", overwrite = TRUE)
writeRaster(rast_temp_proc, "data/rast_temp_proc.tif", overwrite = TRUE)

#### looking at the histograms of the variables ################################

rast_list <- list(
  elevation  = rast_dem_proc,
  ndvi       = rast_ndvi_proc,
  ndwi       = rast_ndwi_proc,
  snowfree   = rast_snowfree_proc,
  slope      = rast_slope_proc,
  hli        = rast_hli_proc,
  tpi        = rast_tpi_proc,
  temp       = rast_temp_proc
)

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

for (v in names(rast_list)) {
  plots <- mp_abiotic[[v]]
  plots <- plots[!is.na(plots)]
  
  land <- values(rast_list[[v]], na.rm = TRUE)
  land <- sample(land, min(5000, length(land)))
  
  dl <- density(land)
  dp <- density(plots)
  
  plot(dl, col = "grey50", lwd = 2, main = v, xlab = v,
       xlim = range(c(dl$x, dp$x)),
       ylim = c(0, max(dl$y, dp$y)))
  lines(dp, col = "steelblue", lwd = 2, lty = 2)
}

par(mfrow = c(1, 1))

# looking at the actual values of the range match
sapply(names(rast_list), function(v) {
  land <- values(rast_list[[v]], na.rm = TRUE)
  p <- mp_abiotic[[v]]; p <- p[!is.na(p)]
  c(land_min = min(land), plot_min = min(p),
    land_max = max(land), plot_max = max(p))
}) |> t() |> round(2)
