#### packages ##################################################################
# install once only - comment out after first run
# devtools::install_github('cjcarlson/embarcadero')

library(embarcadero)
library(pROC)
library(dbarts)
library(raster)
library(parallel)
library(doParallel)
library(foreach)
library(tidyverse)
library(terra)
library(sf)


#### TO DO #####################################################################

# [] Beta support med faktiske værdier og ikke kun threshold

#### setup #####################################################################

set.seed(42)

pred_names <- c("elevation", "slope", "aspect_sin", "aspect_cos", 
                "twi", "temp_predicted", "snowfree")

#### prepare plot data #########################################################

# Get modelable species
modelable_species <- species_frequency |>
  filter(n_plots >= 10) |>
  pull(species)

# Impute missing soil temperature from interpolated raster
imputed_temp <- terra::extract(temp_rast_masked, 
                               plots_sf)[, 2]

# Build predictor + PA matrix
pa_matrix <- species_matrix |>
  dplyr::select(plot_name, all_of(modelable_species)) |>
  mutate(across(-plot_name, ~ as.integer(. > 0))) |>
  left_join(
    abiotic_plot |>
      dplyr::select(plot_name, elevation, slope, aspect_sin, aspect_cos,
                    twi, soil_tem_ave, snowfree) |>
      mutate(temp_predicted = ifelse(is.na(soil_tem_ave),
                                     imputed_temp,
                                     soil_tem_ave)) |>
      dplyr::select(-soil_tem_ave),
    by = "plot_name"
  )

# Verify no NAs in predictors
pa_matrix |>
  dplyr::select(all_of(pred_names)) |>
  summarise(across(everything(), ~sum(is.na(.x)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_na") |>
  filter(n_na > 0)

# Predictor matrix for BART
x_train <- pa_matrix |>
  dplyr::select(all_of(pred_names)) |>
  as.data.frame()

#### read all rasters ##########################################################
dem_rast         <- rast("data/dem_crop.tif")
slope_rast       <- rast("data/slope_crop.tif")
aspect_sin_rast  <- rast("data/aspect_sin_crop.tif")
aspect_cos_rast  <- rast("data/aspect_cos_crop.tif")
twi_rast         <- rast("data/twi_calculated.tif")
ndvi_rast        <- rast("data/ndvi_crop.tif")
temp_rast_masked <- rast("data/temp_predicted_masked.tif")
snowfree_rast    <- rast("data/snow_free_days.tif") |>
  project("EPSG:32622") |>
  crop(aoi) |>
  resample(ndvi_rast)

#### prepare raster stack ######################################################

# Resample all rasters to match ndvi_rast (10m reference)
elev_resamp       <- resample(dem_rast, ndvi_rast)
slope_resamp      <- resample(slope_rast, ndvi_rast)
aspect_sin_resamp <- resample(aspect_sin_rast, ndvi_rast)
aspect_cos_resamp <- resample(aspect_cos_rast, ndvi_rast)
twi_resamp        <- resample(twi_rast, ndvi_rast)
temp_resamp       <- resample(temp_rast_masked, ndvi_rast)

# Stack and name to match predictor names
pred_rast_stack <- c(elev_resamp, slope_resamp, aspect_sin_resamp,
                     aspect_cos_resamp, twi_resamp, temp_resamp,
                     snowfree_rast)

names(pred_rast_stack) <- pred_names

# Convert to old raster format for embarcadero/dbarts
pred_rast_stack_r <- raster::stack(pred_rast_stack)

# Dataframe for prediction — keep NAs, track complete cases
pred_df_r     <- as.data.frame(pred_rast_stack_r, na.rm = FALSE)
complete_idx  <- complete.cases(pred_df_r)

#### BART fitting and spatial projection - parallel ############################

cl <- makeCluster(4)
registerDoParallel(cl)
clusterSetRNGStream(cl, 42)
clusterExport(cl, c("pa_matrix", "x_train", "pred_df_r",
                    "complete_idx", "pred_rast_stack_r", "modelable_species"))

foreach(sp = modelable_species,
        .packages = c("dbarts", "raster"),
        .errorhandling = "pass") %dopar% {
          
          train_data <- cbind(x_train, y = pa_matrix[[sp]]) |> as.data.frame()
          
          model <- dbarts::bart2(
            y ~ .,
            data = train_data,
            keepTrees = TRUE,
            seed = which(modelable_species == sp)
          )
          
          # Predict only on complete cases
          pred_vals <- dbarts:::predict.bart(model, 
                                             newdata = pred_df_r[complete_idx, ])
          prob_vals <- colMeans(pred_vals)
          
          # Insert predictions into full vector with NAs for masked pixels
          prob_full <- rep(NA_real_, nrow(pred_df_r))
          prob_full[complete_idx] <- prob_vals
          
          sp_rast_r <- raster::raster(pred_rast_stack_r[[1]])
          raster::values(sp_rast_r) <- prob_full
          
          filename <- paste0("data/sdm_", gsub(" ", "_", sp), ".tif")
          raster::writeRaster(sp_rast_r, filename, overwrite = TRUE)
          
          filename
        }

stopCluster(cl)

# Load saved rasters back into list
species_rasts <- lapply(modelable_species, function(sp) {
  rast(paste0("data/sdm_", gsub(" ", "_", sp), ".tif"))
})
names(species_rasts) <- modelable_species

species_rasts <- lapply(modelable_species, function(sp) {
  rast(paste0("data/sdm_", gsub(" ", "_", sp), ".tif"))
})
names(species_rasts) <- modelable_species

# Verify they loaded correctly
plot(species_rasts[["Betula nana"]])

#### diagnostics ###############################################################

# Refit models in main session for diagnostics (no parallel - need objects)
bart_models <- list()

for(sp in modelable_species) {
  cat("Fitting BART for diagnostics:", sp, "\n")
  
  train_data <- cbind(x_train, y = pa_matrix[[sp]]) |> as.data.frame()
  
  bart_models[[sp]] <- dbarts::bart2(
    y ~ .,
    data = train_data,
    keepTrees = TRUE,
    seed = which(modelable_species == sp)
  )
}

# Training AUC
auc_results <- map_dfr(modelable_species, function(sp) {
  predicted <- colMeans(
    dbarts:::predict.bart(bart_models[[sp]], 
                          newdata = as.data.frame(x_train))
  )
  observed <- pa_matrix[[sp]]
  
  auc_val <- auc(roc(observed, predicted, quiet = TRUE))
  
  tibble(
    species = sp,
    n_plots = sum(observed),
    auc = round(as.numeric(auc_val), 3)
  )
}) |>
  arrange(desc(auc))
#### explanary ###############################################################

#### variable importance #######################################################

# Refit models in main session for diagnostics
bart_models <- list()

for(sp in modelable_species) {
  cat("Fitting BART for diagnostics:", sp, "\n")
  
  train_data <- cbind(x_train, y = pa_matrix[[sp]]) |> as.data.frame()
  
  bart_models[[sp]] <- dbarts::bart2(
    y ~ .,
    data = train_data,
    keepTrees = TRUE,
    seed = which(modelable_species == sp)
  )
}

# Extract variable importance for all species
varimp_long <- map_dfr(modelable_species, function(sp) {
  vi <- varimp(bart_models[[sp]])
  vi |> mutate(species = sp)
})

# Plot heatmap
ggplot(varimp_long, aes(x = names, y = species, fill = varimps)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Importance") +
  labs(x = "Variable", y = "Species",
       title = "Variable importance across BART SDMs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))


#### variable importance #######################################################
r2_results <- map_dfr(modelable_species, function(sp) {
  predicted <- colMeans(
    dbarts:::predict.bart(bart_models[[sp]], 
                          newdata = as.data.frame(x_train))
  )
  observed <- pa_matrix[[sp]]
  
  # R-squared
  ss_res <- sum((observed - predicted)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  r2 <- 1 - (ss_res / ss_tot)
  
  tibble(
    species = sp,
    r2 = round(r2, 3),
    unexplained = round(1 - r2, 3)
  )
}) |>
  arrange(desc(r2))

# Plot stacked bar of explained vs unexplained
r2_results |>
  pivot_longer(cols = c(r2, unexplained),
               names_to = "component",
               values_to = "value") |>
  mutate(component = factor(component, 
                            levels = c("unexplained", "r2"),
                            labels = c("Unexplained", "Explained")),
         species = factor(species, levels = r2_results$species)) |>
  ggplot(aes(x = value, y = species, fill = component)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Explained" = "#2E75B6", 
                               "Unexplained" = "#D9D9D9")) +
  labs(x = "Proportion of variance", y = NULL,
       fill = NULL,
       title = "Explained vs unexplained variance per species BART SDM") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))
