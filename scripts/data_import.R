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
devtools::install_github('cjcarlson/embarcadero') #sdm also for interpolation
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

#### loading the data ##########################################################

tms <- readRDS("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/tms_pivot.rds") |> 
  clean_names() |> 
  filter(level == "t1_below6cm") |>          # soil temperature only
  filter(month(date) %in% c(6, 7, 8)) |>  # June-August         # June-August
  group_by(plot) |> 
  reframe(
    temp_mean_tms = mean(temp, na.rm = TRUE),
    mean_soilmoisture_tms = mean(moisture_data_raw, na.rm = TRUE)
  ) |> 
  mutate(
    plot_name = toupper(plot),
    vwc_tms = (mean_soilmoisture_tms - min(mean_soilmoisture_tms)) / 
      (max(mean_soilmoisture_tms) - min(mean_soilmoisture_tms)) * 100
  )

summary(tms)

BioBasis_Nuuk_PhenologyPlots_Microclimate_2025 <- read_excel("~/Library/CloudStorage/OneDrive-GrønlandsNaturinstitut/General - BioBasis/03_GEM_Database/Datafiler excel/BioBasis_Nuuk_PhenologyPlots_Microclimate_2025.xlsx")
BioBasis_Nuuk_CFlux_Microclimate_2025 <- read_excel("~/Library/CloudStorage/OneDrive-GrønlandsNaturinstitut/General - BioBasis/03_GEM_Database/Datafiler excel/BioBasis_Nuuk_CFlux_Microclimate_2025.xlsx")

tms_biobasis <- BioBasis_Nuuk_CFlux_Microclimate_2025 |> 
  bind_rows(BioBasis_Nuuk_PhenologyPlots_Microclimate_2025) |> 
  filter(month(Date) %in% c(6, 7, 8)) |>  # June-August
  group_by(Plot, Latitude, Longitude) |>
  reframe(
    temp_mean_tms = mean(Temp_6cmbel, na.rm = TRUE),
    mean_soilmoisture_tms = mean(Raw_soil_moisture, na.rm = TRUE)
  ) |>
  mutate(
    vwc_tms = (mean_soilmoisture_tms - min(mean_soilmoisture_tms)) / 
      (max(mean_soilmoisture_tms) - min(mean_soilmoisture_tms)) * 100
  )

summary(tms_biobasis)

samples_qgis <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/samples_qgis.csv") |> 
  dplyr::select(plot, X,Y,elevation, ndvi, ndwi) |> 
  mutate(plot_name = plot) |> 
  left_join(tms, by = "plot_name")

names(samples_qgis)

df_raw <- read_csv("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/r_data/future_changes_data/data/samples.csv", col_types = cols(Date = col_datetime(format = "%m/%d/%Y %H.%M"))) |> 
  clean_names() |> 
  mutate(rowid = row_number(), plot_name = toupper(plot_name))

df_cover <- df_raw |> 
 mutate(across(ends_with("_bb"), bb_to_cover)) |> 
  mutate(total_cover = rowSums(across(ends_with("_bb")), na.rm = TRUE))

summary(df_cover)

#### cobining all tms ##########################################################

# Bind and normalise vwc together
tms_combined <- bind_rows(tms_own, tms_biobasis) |>
  mutate(
    vwc_tms = (mean_soilmoisture_tms - min(mean_soilmoisture_tms)) / 
      (max(mean_soilmoisture_tms) - min(mean_soilmoisture_tms)) * 100
  )

summary(tms_combined)

nrow(tms_combined)

tms_combined_sf <- tms_combined |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) |>
  st_transform(32622)

plot(ndvi_rast)
plot(st_geometry(tms_combined_sf), add = TRUE, pch = 16, col = "red", cex = 0.8)

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
  dplyr::select(plot_name, species_name, cover) |>
  pivot_wider(
    names_from = species_name,
    values_from = cover,
    values_fill = 0
  )

sp_cols <- species_matrix |> dplyr::select(-plot_name)

abiotic_plot <- df_cover |>
  left_join(species_matrix |> dplyr::select(plot_name), by = "plot_name") |>
  mutate(
    richness = rowSums(sp_cols > 0),
    shannon = vegan::diversity(sp_cols, index = "shannon")
  )

#### final abiotic df###########################################################

abiotic_plot <- abiotic_plot |> 
  dplyr::select(plot_name, veg_height_ave, bare_ground_bb, x, y, total_cover, richness, shannon, soil_moi_ave, soil_tem_ave) |> 
  left_join(tms, by = "plot_name") |> 
  dplyr::select(-plot)
  
summary(abiotic_plot)

# Find the plot with NA temp
na_plot <- abiotic_plot |> 
  filter(is.na(soil_tem_ave))

# Convert to sf and extract from raster
na_plot_sf <- na_plot |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(32622)

# Extract interpolated value
imputed_temp <- terra::extract(temp_rast, na_plot_sf)[, 2]
imputed_temp

# Fill NA in abiotic_plot
abiotic_plot <- abiotic_plot |>
  mutate(temp_predicted = ifelse(is.na(soil_tem_ave), 
                                 imputed_temp, 
                                 soil_tem_ave))

# Verify no more NAs
sum(is.na(abiotic_plot$temp_predicted))

# Extract TMS logger plots with coordinates from abiotic_plot
tms_own <-  abiotic_plot |>
  filter(!is.na(temp_mean_tms)) |>
  dplyr::select(plot_name, x, y, temp_mean_tms, mean_soilmoisture_tms) |>
  rename(Longitude = x, Latitude = y)

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

plots_sf <- abiotic_plot |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(32622)

aoi <- plots_sf |>
  st_bbox() |>
  st_as_sfc() |>
  st_buffer(50) |>
  vect()  # convert to terra format for cropping



#### importing ndvi ############################################################

ndvi_rast <- rast("data/ndvi_export_2025.tif") |> 
  crop(aoi)

plot(ndvi_rast)
summary(ndvi_rast)
print(ndvi_rast)

abiotic_plot <- abiotic_plot |>
#  select(-"ndvi") |>
  mutate(ndvi = terra::extract(ndvi_rast, plots_sf)[, 2])

abiotic_plot |> 
  dplyr::select(plot_name, ndvi) |> 
  summary()

#### importing elevation ############################################################

dem_rast <- rast("data/elevation_arcticdem-30_32622.tif") |> 
  crop(aoi)

plot(dem_rast)
summary(dem_rast)
print(dem_rast)

abiotic_plot <- abiotic_plot |>
  mutate(elevation = terra::extract(dem_rast, plots_sf)[, 2])

abiotic_plot |> 
  dplyr::select(plot_name, elevation) |> 
  summary()

#### slope #####################################################################
slope_rast <- terrain(dem_rast, v = "slope", unit = "degrees")

abiotic_plot <- abiotic_plot |>
  mutate(slope = terra::extract(slope_rast, plots_sf)[, 2])

abiotic_plot |> 
  dplyr::select(plot_name, slope) |> 
  summary()

#### aspect ####################################################################

aspect_rast <- terrain(dem_rast, v = "aspect", unit = "degrees") |> 
  crop(aoi)

abiotic_plot <- abiotic_plot |>
  mutate(
    aspect_raw = terra::extract(aspect_rast, plots_sf)[, 2],
    aspect_sin = sin(aspect_raw * pi / 180),
    aspect_cos = cos(aspect_raw * pi / 180)
  )

abiotic_plot |> 
  dplyr::select(plot_name, aspect_raw, aspect_sin, aspect_cos) |> 
  summary()

aspect_cos_rast <- cos(aspect_rast * pi / 180) |> 
  crop(aoi)
aspect_sin_rast <- sin(aspect_rast * pi / 180) |> 
  crop(aoi)
#### calculating twi #########################################################

writeRaster(dem_rast, "data/dem_crop.tif", overwrite = TRUE)

wbt_fill_depressions("data/dem_crop.tif", "data/dem_filled.tif")
wbt_d8_flow_accumulation("data/dem_filled.tif", "data/sca.tif")
wbt_slope("data/dem_filled.tif", "data/slope_wb.tif", units = "degrees")
wbt_wetness_index(
  sca = "data/sca.tif",
  slope = "data/slope_wb.tif",
  output = "data/twi_calculated.tif"
)

twi_rast <- rast("data/twi_calculated.tif")

abiotic_plot <- abiotic_plot |>
  mutate(twi = terra::extract(twi_rast, plots_sf)[, 2])

abiotic_plot |> 
  dplyr::select(plot_name, twi) |> 
  summary()


#### checking the nas ##########################################################

abiotic_plot |> 
  filter(is.na(elevation) | is.na(slope) | is.na(aspect_raw)) |> 
  dplyr::select(plot_name, x, y, elevation, slope, aspect_raw, twi, ndvi)

summary(abiotic_plot)

#### writing all data files ####################################################
saveRDS(abiotic_plot, "data/abiotic_plot.rds")
saveRDS(species_matrix, "data/species_matrix.rds")
saveRDS(species_long, "data/species_long.rds")
writeRaster(dem_rast, "data/dem_crop.tif", overwrite = TRUE)
writeRaster(twi_rast, "data/twi_calculated.tif", overwrite = TRUE)
writeRaster(ndvi_rast, "data/ndvi_crop.tif", overwrite = TRUE)
writeRaster(slope_rast, "data/slope_crop.tif", overwrite = TRUE)
writeRaster(aspect_rast, "data/aspect_crop.tif", overwrite = TRUE)

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

##### temperature interpolation ###########################################

# Step 1: fit linear model on combined logger data
temp_lm <- lm(temp_mean_tms ~ ndvi + aspect_cos + aspect_sin + slope, 
              data = tms_combined)

temp_lm2 <- lm(temp_mean_tms ~ ndvi + aspect_cos + aspect_sin, 
              data = tms_combined)

summary(temp_lm)

# Add residuals to combined logger data
 tms_combined <- tms_combined |>
  mutate(temp_resid = residuals(temp_lm2))

# Convert to sf for variogram
tms_combined_sf <- tms_combined |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) |>
  st_transform(32622)

# Compute variogram of residuals
vgm_temp <- variogram(temp_resid ~ 1, data = tms_combined_sf)
plot(vgm_temp)

# Resample aspect layers to match ndvi resolution and extent
aspect_cos_rast <- resample(aspect_cos_rast, ndvi_rast)
aspect_sin_rast <- resample(aspect_sin_rast, ndvi_rast)

# Now stack
pred_stack <- c(ndvi_rast, aspect_cos_rast, aspect_sin_rast)
names(pred_stack) <- c("ndvi", "aspect_cos", "aspect_sin")

# Project
temp_rast <- predict(pred_stack, temp_lm2)
plot(temp_rast)

writeRaster(temp_rast, "data/temp_predicted_rast.tif", overwrite = TRUE)

temp_rast_masked <- mask(temp_rast, ndvi_rast < 0.1, maskvalue = TRUE)

plot(temp_rast_masked)
writeRaster(temp_rast_masked, "data/temp_predicted_rast.tif", overwrite = TRUE)

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
    dplyr::select(ndvi, twi, elevation, slope, aspect_sin, aspect_cos) |>
    as.data.frame(),
  y.train = abiotic_plot$soil_moi_ave,
  keeptrees = TRUE
)

summary(moist_bart)

#### BART ####

predictors <- c("elevation", "slope", "aspect_sin", "aspect_cos", 
                "twi", "ndvi", "soil_tem_ave")

plot_predictors <- abiotic_plot |>
  dplyr::select(all_of(predictors), plot_name)

modelable_species <- species_frequency |>
  filter(n_plots >= 10) |>
  pull(species)

modelable_species

# Create binary presence/absence for modelable species only
pa_matrix <- species_matrix |>
  dplyr::select(plot_name, all_of(modelable_species)) |>
  mutate(across(-plot_name, ~ as.integer(. > 0))) |>
  left_join(plot_predictors, by = "plot_name")

glimpse(pa_matrix)

pa_matrix <- pa_matrix |>
  mutate(temp_predicted = terra::extract(temp_rast_masked, plots_sf)[, 2]) |>
  dplyr::select(-soil_tem_ave)

pa_matrix |>
  dplyr::select(elevation, slope, aspect_sin, aspect_cos, twi, ndvi, temp_predicted) |>
  summarise(across(everything(), ~sum(is.na(.x)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_na") |>
  filter(n_na > 0)

### BART loop #################

# Define predictor names
pred_names <- c("elevation", "slope", "aspect_sin", "aspect_cos", 
                "twi", "ndvi", "temp_predicted")

# Prepare predictor dataframe
x_train <- pa_matrix |>
  dplyr::select(all_of(pred_names)) |>
  as.data.frame()

# Run BART SDM for each species and store results
bart_models <- list()

for(sp in modelable_species) {
  cat("Fitting BART for:", sp, "\n")
  
  y_train <- pa_matrix[[sp]]
  
  bart_models[[sp]] <- bart(
    x.train = x_train,
    y.train = y_train,
    keeptrees = TRUE
  )
}

names(bart_models)
length(bart_models)

#checking what specis is missing

modelable_species[!modelable_species %in% names(bart_models)]
##### diagnostics ####

auc_results <- map_dfr(modelable_species, function(sp) {
  predicted <- pnorm(colMeans(bart_models[[sp]]$yhat.train))
  observed <- pa_matrix[[sp]]
  
  auc_val <- auc(roc(observed, predicted, quiet = TRUE))
  
  tibble(
    species = sp,
    n_plots = sum(observed),
    auc = round(as.numeric(auc_val), 3)
  )
}) |>
  arrange(desc(auc))

print(auc_results, n = Inf)

##### 5 fold validation ########################################################

cv_results <- map_dfr(modelable_species, function(sp) {
  y <- pa_matrix[[sp]]
  
  # Stratified folds - ensure presences and absences in each fold
  pres_idx <- which(y == 1)
  abs_idx <- which(y == 0)
  
  pres_folds <- sample(rep(1:5, length.out = length(pres_idx)))
  abs_folds <- sample(rep(1:5, length.out = length(abs_idx)))
  
  folds <- numeric(100)
  folds[pres_idx] <- pres_folds
  folds[abs_idx] <- abs_folds
  
  fold_auc <- map_dbl(1:5, function(k) {
    train_idx <- folds != k
    test_idx <- folds == k
    
    cv_model <- bart(
      x.train = x_train[train_idx, ],
      y.train = y[train_idx],
      x.test = x_train[test_idx, ],
      keeptrees = FALSE
    )
    
    pred <- pnorm(colMeans(cv_model$yhat.test))
    auc(roc(y[test_idx], pred, quiet = TRUE)) |> as.numeric()
  })
  
  tibble(species = sp, cv_auc = mean(fold_auc), sd_auc = sd(fold_auc))
}) |>
  arrange(desc(cv_auc))

print(cv_results, n = Inf)

cv_results <- cv_results |>
  mutate(assessment = case_when(
    cv_auc >= 0.90 ~ "Excellent",
    cv_auc >= 0.80 ~ "Good",
    cv_auc >= 0.70 ~ "Acceptable",
    cv_auc >= 0.60 ~ "Poor",
    TRUE           ~ "Fail"
  ))

print(cv_results, n = Inf)

#### spacially explicit ###

# Resample all rasters to match ndvi_rast (10m reference)
elev_resamp <- resample(dem_rast, ndvi_rast)
slope_resamp <- resample(slope_rast, ndvi_rast)
aspect_sin_resamp <- resample(aspect_sin_rast, ndvi_rast)
aspect_cos_resamp <- resample(aspect_cos_rast, ndvi_rast)
twi_resamp <- resample(twi_rast, ndvi_rast)
temp_resamp <- resample(temp_rast_masked, ndvi_rast)

# Stack all predictors
pred_rast_stack <- c(elev_resamp, slope_resamp, aspect_sin_resamp, 
                     aspect_cos_resamp, twi_resamp, ndvi_rast, temp_resamp)

# Name layers to match predictor names in BART models
names(pred_rast_stack) <- c("elevation", "slope", "aspect_sin", "aspect_cos",
                            "twi", "ndvi", "temp_predicted")

print(pred_rast_stack)
plot(pred_rast_stack)

# Convert raster stack to dataframe for prediction
pred_df <- as.data.frame(pred_rast_stack, xy = TRUE, na.rm = FALSE)

# Get just the predictor columns
pred_only <- pred_df |>
  dplyr::select(all_of(pred_names)) |>
  as.data.frame()


pred_rast_stack_r <- raster::stack(pred_rast_stack)

# Test prediction
test_pred <- embarcadero:::predict2.bart(bart_models[["Coptis trifolia"]], 
                                         x.layers = pred_rast_stack_r)
dim(test_pred)

test_rast <- rast(test_pred[,,1], 
                  extent = ext(ndvi_rast),
                  crs = crs(ndvi_rast))

plot(test_rast)
summary(test_rast)

##### testing loop

species_rasts <- list()

for(sp in modelable_species) {
  cat("Projecting:", sp, "\n")
  
  pred_vals <- embarcadero:::predict2.bart(bart_models[[sp]], 
                                           x.layers = pred_rast_stack_r)
  
  # Convert to terra rast using first layer as template
  sp_rast_r <- raster::raster(pred_rast_stack_r[[1]])
  raster::values(sp_rast_r) <- pred_vals[[1]][]
  
  species_rasts[[sp]] <- rast(sp_rast_r)
}

# Save all species rasters
for(sp in modelable_species) {
  filename <- paste0("data/sdm_", gsub(" ", "_", sp), ".tif")
  writeRaster(species_rasts[[sp]], filename, overwrite = TRUE)
}

#### diagnostics ####################################################

cv_results <- cv_auc

diagnostics <- map_dfr(modelable_species, function(sp) {
  
  # Training predictions
  predicted_train <- colMeans(pnorm(bart_models[[sp]]$yhat.train))
  observed <- pa_matrix[[sp]]
  
  # Training AUC
  train_auc <- auc(roc(observed, predicted_train, quiet = TRUE)) |> as.numeric()
  
  # Prevalence
  prevalence <- mean(observed)
  
  # TSS at optimal threshold
  roc_obj <- roc(observed, predicted_train, quiet = TRUE)
  best_thresh <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
  tss <- best_thresh$sensitivity + best_thresh$specificity - 1
  
  tibble(
    species = sp,
    n_presences = sum(observed),
    prevalence = round(prevalence, 2),
    train_auc = round(train_auc, 3),
    cv_auc = round(cv_results$cv_auc[cv_results$species == sp], 3),
    cv_auc_sd = round(cv_results$sd_auc[cv_results$species == sp], 3),
    tss = round(tss, 3)
  )
}) |>
  arrange(desc(cv_auc))

print(diagnostics, n = Inf)

# TO DO 
# check that AUC is on the k-fold validation data
# check number of species where that are good

learning_curve <- map_dfr(modelable_species, function(sp) {
  y <- pa_matrix[[sp]]
  
  sample_sizes <- c(20, 30, 40, 50, 60, 70, 80, 90, 100)
  
  map_dfr(sample_sizes, function(n) {
    aucs <- map(1:5, function(i) {
      tryCatch({
        # Stratified sampling to ensure both classes in train and test
        pres_idx <- which(y == 1)
        abs_idx <- which(y == 0)
        
        n_pres <- round(n * mean(y))
        n_abs <- n - n_pres
        
        train_idx <- c(
          sample(pres_idx, min(n_pres, length(pres_idx))),
          sample(abs_idx, min(n_abs, length(abs_idx)))
        )
        test_idx <- setdiff(1:100, train_idx)
        
        # Skip if test set has only one class
        if(length(unique(y[test_idx])) < 2) return(NULL)
        
        m <- bart(
          x.train = x_train[train_idx, ],
          y.train = y[train_idx],
          x.test = x_train[test_idx, ],
          keeptrees = FALSE
        )
        
        pred <- pnorm(colMeans(m$yhat.test))
        auc(roc(y[test_idx], pred, quiet = TRUE)) |> as.numeric()
      }, error = function(e) NULL)
    }) |> 
      compact() |>  # remove NULLs
      unlist()
    
    tibble(species = sp, n_train = n, mean_auc = mean(aucs))
  })
})

# Plot
ggplot(learning_curve, aes(x = n_train, y = mean_auc, colour = species)) +
  geom_line() +
  geom_vline(xintercept = 100, linetype = "dashed") +
  geom_vline(xintercept = 120, linetype = "dotted", colour = "red") +
  labs(x = "Training sample size", y = "Mean AUC",
       title = "Learning curves") +
  theme_minimal()

varimp_results <- map_dfr(modelable_species, function(sp) {
  varimp(bart_models[[sp]]) |>
    mutate(species = sp)
}) |>
  dplyr::select(species, names, varimps) |>
  pivot_wider(names_from = names, values_from = varimps) |>
  arrange(species)

print(varimp_results, n = Inf)

varimp_long <- map_dfr(modelable_species, function(sp) {
  varimp(bart_models[[sp]]) |>
    mutate(species = sp)
})

ggplot(varimp_long, aes(x = names, y = species, fill = varimps)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "Variable", y = "Species", fill = "Importance",
       title = "Variable importance across BART SDMs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
