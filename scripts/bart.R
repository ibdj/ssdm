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

#### setup #####################################################################

set.seed(42)

pred_names <- c("elevation", "slope", "aspect_sin", "aspect_cos", 
                "twi", "temp_predicted")

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
                    twi, ndvi, soil_tem_ave) |>
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

#### prepare raster stack ######################################################

# Resample all rasters to match ndvi_rast (10m reference)
elev_resamp        <- resample(dem_rast, ndvi_rast)
slope_resamp       <- resample(slope_rast, ndvi_rast)
aspect_sin_resamp  <- resample(aspect_sin_rast, ndvi_rast)
aspect_cos_resamp  <- resample(aspect_cos_rast, ndvi_rast)
twi_resamp         <- resample(twi_rast, ndvi_rast)
temp_resamp        <- resample(temp_rast_masked, ndvi_rast)

# Stack and name to match predictor names
pred_rast_stack <- c(elev_resamp, slope_resamp, aspect_sin_resamp,
                     aspect_cos_resamp, twi_resamp, temp_resamp)

names(pred_rast_stack) <- pred_names

# Convert to old raster format for embarcadero/dbarts
pred_rast_stack_r <- raster::stack(pred_rast_stack)

# Dataframe for prediction — keep NAs, track complete cases
pred_df_r <- as.data.frame(pred_rast_stack_r, na.rm = FALSE)
complete_idx <- complete.cases(pred_df_r)

#### BART fitting and spatial projection - parallel ############################

cl <- makeCluster(4)
registerDoParallel(cl)
clusterSetRNGStream(cl, 42)
clusterExport(cl, c("pa_matrix", "x_train", "pred_df_r",
                    "complete_idx", "pred_rast_stack_r", "modelable_species"))

foreach(sp = modelable_species,
        .packages = c("dbarts", "raster", "terra"),
        .errorhandling = "pass") %dopar% {
          
          train_data <- cbind(x_train, y = pa_matrix[[sp]]) |> as.data.frame()
          model <- dbarts::bart2(y ~ .,
                                 data = train_data,
                                 keepTrees = TRUE,
                                 seed = which(modelable_species == sp))
          
          pred_vals <- dbarts:::predict.bart(model, newdata = pred_df_r)
          prob_vals <- colMeans(pred_vals)
          
          prob_full <- rep(NA_real_, nrow(pred_df_r))
          prob_full[complete_idx] <- prob_vals
          
          sp_rast_r <- raster::raster(pred_rast_stack_r[[1]])
          raster::values(sp_rast_r) <- prob_full
          
          # Save directly from worker - don't return raster object
          filename <- paste0("data/sdm_", gsub(" ", "_", sp), ".tif")
          raster::writeRaster(sp_rast_r, filename, overwrite = TRUE)
          
          filename  # return path as confirmation
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

auc_results <- auc_results |>
  mutate(assessment = case_when(
    auc >= 0.90 ~ "Excellent",
    auc >= 0.80 ~ "Good",
    auc >= 0.70 ~ "Acceptable",
    auc >= 0.60 ~ "Poor",
    TRUE        ~ "Fail"
  ))

print(auc_results, n = Inf)

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

print(cv_results, n = Inf)

#### spacially explicit ########################################################
################################################################################

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

##### testing loop #############################################################
################################################################################

cl <- makeCluster(4)
registerDoParallel(cl)

# Export required objects to workers
clusterExport(cl, c("bart_models", "pred_rast_stack_r", "modelable_species"))

species_rasts_list <- foreach(sp = modelable_species,
                              .packages = c("embarcadero", "raster", "terra"),
                              .errorhandling = "pass") %dopar% {
                                
                                pred_vals <- embarcadero:::predict2.bart(bart_models[[sp]], 
                                                                         x.layers = pred_rast_stack_r)
                                
                                sp_rast_r <- raster::raster(pred_rast_stack_r[[1]])
                                raster::values(sp_rast_r) <- raster::values(pred_vals[[1]])
                                
                                rast(sp_rast_r)
                              }

stopCluster(cl)

# Restore named list
names(species_rasts_list) <- modelable_species
species_rasts <- species_rasts_list

# Save all species rasters
for(sp in modelable_species) {
  filename <- paste0("data/sdm_", gsub(" ", "_", sp), ".tif")
  writeRaster(species_rasts[[sp]], filename, overwrite = TRUE)
}

#### diagnostics ###############################################################
################################################################################

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
