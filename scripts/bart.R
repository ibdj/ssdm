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



#### setup #####################################################################

set.seed(42)

pred_names <- c("elevation", "slope", "aspect_sin", "aspect_cos", 
                "ndwi", "temp_predicted", "snowfree")

#### prepare plot data #########################################################

read_rds("data/species_frequency.rds")

# Get modelable species
modelable_species <- species_frequency |>
  filter(n_plots >= 10) |>
  pull(species)

temp_rast_masked <- rast("data/temp_predicted_masked.tif")

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
aspect_sin_rast  <- rast("data/aspect_cos_rast.tif")
aspect_cos_rast  <- rast("data/aspect_sin_rast.tif")
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

##### 5 fold validation ########################################################

cl <- makeCluster(4)
registerDoParallel(cl)
clusterSetRNGStream(cl, 42)
clusterExport(cl, c("pa_matrix", "x_train", "modelable_species"))

cv_results <- foreach(sp = modelable_species,
                      .packages = c("dbarts", "pROC"),
                      .combine = rbind,
                      .errorhandling = "pass") %dopar% {
                        
                        y <- pa_matrix[[sp]]
                        
                        # Stratified folds
                        pres_idx <- which(y == 1)
                        abs_idx  <- which(y == 0)
                        folds <- numeric(100)
                        folds[pres_idx] <- sample(rep(1:5, length.out = length(pres_idx)))
                        folds[abs_idx]  <- sample(rep(1:5, length.out = length(abs_idx)))
                        
                        fold_auc <- sapply(1:5, function(k) {
                          train_idx <- which(folds != k)
                          test_idx  <- which(folds == k)
                          
                          train_data <- cbind(x_train[train_idx, ], y = y[train_idx]) |> 
                            as.data.frame()
                          
                          cv_model <- dbarts::bart2(
                            y ~ .,
                            data = train_data,
                            keepTrees = TRUE,
                            seed = which(modelable_species == sp) * k
                          )
                          
                          pred <- colMeans(
                            dbarts:::predict.bart(cv_model, 
                                                  newdata = as.data.frame(x_train[test_idx, ]))
                          )
                          
                          as.numeric(auc(roc(y[test_idx], pred, quiet = TRUE)))
                        })
                        
                        data.frame(species = sp, cv_auc = mean(fold_auc), sd_auc = sd(fold_auc))
                      }

stopCluster(cl)

cv_results <- cv_results |>
  arrange(desc(cv_auc)) |>
  mutate(assessment = case_when(
    cv_auc >= 0.90 ~ "Excellent",
    cv_auc >= 0.80 ~ "Good",
    cv_auc >= 0.70 ~ "Acceptable",
    cv_auc >= 0.60 ~ "Poor",
    TRUE           ~ "Fail"
  ))

cv_results |> 
  tibble::as_tibble() |> 
  print(n = Inf)

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

varimp_long <- map_dfr(modelable_species, function(sp) {
  vc <- bart_models[[sp]]$varcount
  
  # Average across chains (dim 1) and iterations (dim 2)
  counts <- apply(vc, 3, mean)
  
  tibble(
    variable = pred_names,
    importance = counts / sum(counts),
    species = sp
  )
})

head(varimp_long)

ggplot(varimp_long, aes(x = variable, y = species, fill = importance)) +
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

### expl. vs. unexpl. ###############

# R2 drop when each variable is excluded
varimp_r2 <- map_dfr(modelable_species, function(sp) {
  observed <- pa_matrix[[sp]]
  
  # Full model R2
  full_pred <- colMeans(dbarts:::predict.bart(bart_models[[sp]], 
                                              newdata = as.data.frame(x_train)))
  ss_tot <- sum((observed - mean(observed))^2)
  r2_full <- 1 - sum((observed - full_pred)^2) / ss_tot
  
  # R2 drop per variable
  map_dfr(pred_names, function(var) {
    x_reduced <- x_train
    x_reduced[[var]] <- mean(x_train[[var]], na.rm = TRUE)  # replace with mean
    
    red_pred <- colMeans(dbarts:::predict.bart(bart_models[[sp]], 
                                               newdata = as.data.frame(x_reduced)))
    r2_reduced <- 1 - sum((observed - red_pred)^2) / ss_tot
    
    tibble(species = sp, variable = var, 
           r2_drop = max(0, r2_full - r2_reduced))
  })
})


# Recalculate unexplained from r2_results
unexplained <- r2_results |>
  mutate(variable = "Unexplained",
         r2_drop = 1 - r2)

# Combine
plot_data <- varimp_r2 |>
  dplyr::select(species, variable, r2_drop) |>
  bind_rows(dplyr::select(unexplained, species, variable, r2_drop)) |>
  mutate(
    species = factor(species, levels = r2_results$species),
    variable = factor(variable, levels = c("Unexplained", rev(pred_names)))
  )

# Plot
ggplot(plot_data, aes(x = r2_drop, y = species, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c(
      "Unexplained"    = "#D9D9D9",
      "elevation"      = "#1f77b4",
      "slope"          = "#ff7f0e",
      "aspect_sin"     = "#2ca02c",
      "aspect_cos"     = "#d62728",
      "twi"            = "#9467bd",
      "temp_predicted" = "#8c564b",
      "snowfree"       = "#e377c2"
    )
  ) +
  labs(x = "Proportion of variance", y = NULL,
       fill = "Variable",
       title = "Explained variance partitioned by predictor") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))


# VERSION TWO ##################################################################

# Normalise r2_drop to sum to r2_full per species
varimp_r2_norm <- varimp_r2 |>
  left_join(r2_results, by = "species") |>
  group_by(species) |>
  mutate(r2_drop_norm = (r2_drop / sum(r2_drop)) * r2) |>
  ungroup()

# Recalculate unexplained
unexplained <- r2_results |>
  mutate(variable = "Unexplained",
         r2_drop_norm = 1 - r2)

# Combine
plot_data <- varimp_r2_norm |>
  dplyr::select(species, variable, r2_drop_norm) |>
  bind_rows(dplyr::select(unexplained, species, variable, r2_drop_norm)) |>
  mutate(
    species = factor(species, levels = r2_results$species),
    variable = factor(variable, levels = c("Unexplained", rev(pred_names)))
  )

# Plot
ggplot(plot_data, aes(x = r2_drop_norm, y = species, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(
    values = c(
      "Unexplained"    = "#D9D9D9",
      "elevation"      = "#1f77b4",
      "slope"          = "#ff7f0e",
      "aspect_sin"     = "#2ca02c",
      "aspect_cos"     = "#d62728",
      "twi"            = "#9467bd",
      "temp_predicted" = "#8c564b",
      "snowfree"       = "#e377c2"
    )
  ) +
  labs(x = "Proportion of variance", y = NULL,
       fill = "Variable",
       title = "Explained variance partitioned by predictor") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))

##### variables ###########

# Pivot x_train to long format for faceted plotting
x_train_long <- x_train |>
  pivot_longer(everything(), 
               names_to = "variable", 
               values_to = "value")

# Histogram of plot distribution across each predictor
ggplot(x_train_long, aes(x = value)) +
  geom_histogram(bins = 20, fill = "#2E75B6", colour = "white") +
  facet_wrap(~ variable, scales = "free") +
  labs(x = "Value", y = "Number of plots",
       title = "Distribution of plots across environmental predictors") +
  theme_minimal()
