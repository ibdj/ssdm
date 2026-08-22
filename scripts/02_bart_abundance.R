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
library(doParallel)

#### reading rds ###############################################################

rs <- readRDS("~/OneDrive - Aarhus universitet/MappingPlants/02 Modelling future changes/ssdm/data/rs.rds")
species_matrix_cover <- readRDS("~/OneDrive - Aarhus universitet/MappingPlants/02 Modelling future changes/ssdm/data/species_matrix_cover.rds")
species_frequency <- readRDS("~/OneDrive - Aarhus universitet/MappingPlants/02 Modelling future changes/ssdm/data/species_frequency.rds")
bart_models_pa <- readRDS("~/OneDrive - Aarhus universitet/MappingPlants/02 Modelling future changes/ssdm/data/bart_models_pa.rds")

#### setup #####################################################################
set.seed(42)

pred_stack <- rast("data/pred_stack.tif")
pred_names <- names(pred_stack)

modelable_species <- readRDS("data/species_frequency.rds") |>
  dplyr::filter(n_plots >= 10) |>
  dplyr::pull(taxon)

cover_matrix <- readRDS("data/species_matrix_cover.rds") |>
  as.data.frame() |>
  tibble::rownames_to_column("plot") |>
  dplyr::select(plot, dplyr::all_of(modelable_species)) |>
  dplyr::left_join(sf::st_drop_geometry(rs), by = "plot")

x_train_all <- cover_matrix |>
  dplyr::select(dplyr::all_of(pred_names)) |>
  as.data.frame()

# prediction df (same as PA script)
pred_rast_stack_r <- raster::stack(pred_stack)
pred_df_r    <- as.data.frame(pred_rast_stack_r, na.rm = FALSE)
complete_idx <- complete.cases(pred_df_r)

#### hurdle fitting and projection - parallel ##################################
cl <- makeCluster(4)
registerDoParallel(cl)
clusterSetRNGStream(cl, 42)
clusterExport(cl, c("cover_matrix", "x_train_all", "pred_df_r",
                    "complete_idx", "pred_rast_stack_r", "modelable_species"))

foreach(sp = modelable_species,
        .packages = c("dbarts", "raster"),
        .errorhandling = "pass") %dopar% {
          
          # conditional data: presence plots only
          pres <- which(cover_matrix[[sp]] > 0)
          train_data <- cbind(x_train_all[pres, ], y = cover_matrix[[sp]][pres]) |>
            as.data.frame()
          
          model <- dbarts::bart2(
            y ~ .,
            data = train_data,
            keepTrees = TRUE,
            seed = which(modelable_species == sp)
          )
          
          invisible(model$fit$state)
          saveRDS(model, paste0("data/bart_cov_", gsub(" ", "_", sp), ".rds"))
          
          # conditional cover over AOI
          cov_vals <- colMeans(dbarts:::predict.bart(model,
                                                     newdata = pred_df_r[complete_idx, ]))
          cov_full <- rep(NA_real_, nrow(pred_df_r))
          cov_full[complete_idx] <- pmax(cov_vals, 0)          # clamp negatives
          
          cov_rast <- raster::raster(pred_rast_stack_r[[1]])
          raster::values(cov_rast) <- cov_full
          
          # hurdle = PA prob × conditional cover
          pa_rast     <- raster::raster(paste0("data/sdm_", gsub(" ", "_", sp), ".tif"))
          hurdle_rast <- pa_rast * cov_rast
          
          raster::writeRaster(cov_rast,
                              paste0("data/cov_", gsub(" ", "_", sp), ".tif"), overwrite = TRUE)
          raster::writeRaster(hurdle_rast,
                              paste0("data/hurdle_", gsub(" ", "_", sp), ".tif"), overwrite = TRUE)
          
          c(sp, length(pres))
        }
stopCluster(cl)


#############################
#### running PA BART first (for hurdle apporach due to 0-inflation) ############

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
plot(trim(species_rasts[["Betula nana"]]))

#### fitting and running cover BART ############################################
# Predictor matrix unchanged from pa data — reuse existing x_train

# Response: cover instead of binary presence
# (in the foreach loop, replace:)
#   train_data <- cbind(x_train, y = pa_matrix[[sp]]) |> as.data.frame()
# with:
train_data <- cbind(x_train, y = cover_matrix[[sp]]) |> as.data.frame()

model <- dbarts::bart2(
  y ~ .,
  data = train_data,
  keepTrees = TRUE,
  seed = which(modelable_species == sp)
)

pred_check <- colMeans(dbarts:::predict.bart(model, newdata = as.data.frame(x_train)))
summary(pred_check)

cl <- makeCluster(4)
registerDoParallel(cl)
clusterSetRNGStream(cl, 42)
clusterExport(cl, c("cover_matrix", "x_train", "pred_df_r",
                    "complete_idx", "pred_rast_stack_r", "modelable_species"))

foreach(sp = modelable_species,
        .packages = c("dbarts", "raster"),
        .errorhandling = "pass") %dopar% {
          
          train_data <- cbind(x_train, y = cover_matrix[[sp]]) |> as.data.frame()
          
          model <- dbarts::bart2(
            y ~ .,
            data = train_data,
            keepTrees = TRUE,
            seed = which(modelable_species == sp)
          )
          
          pred_vals <- dbarts:::predict.bart(model, 
                                             newdata = pred_df_r[complete_idx, ])
          cover_vals <- colMeans(pred_vals)
          
          cover_full <- rep(NA_real_, nrow(pred_df_r))
          cover_full[complete_idx] <- cover_vals
          
          sp_rast_r <- raster::raster(pred_rast_stack_r[[1]])
          raster::values(sp_rast_r) <- cover_full
          sp_rast_r[sp_rast_r < 0] <- 0 #making sure that bart doesnt extrapolate below 0
          
          filename <- paste0("data/sdm_cover_", gsub(" ", "_", sp), ".tif")
          raster::writeRaster(sp_rast_r, filename, overwrite = TRUE)
          
          filename
        }

stopCluster(cl)

species_rasts_cover <- lapply(modelable_species, function(sp) {
  rast(paste0("data/sdm_cover_", gsub(" ", "_", sp), ".tif"))
})
names(species_rasts_cover) <- modelable_species

plot(trim(species_rasts_cover[["Betula nana"]]))

#### 5-fold cross validation with RMSE #########################################

cl <- makeCluster(4)
registerDoParallel(cl)
clusterSetRNGStream(cl, 42)
clusterExport(cl, c("cover_matrix", "x_train", "modelable_species"))

cv_results_cover <- foreach(sp = modelable_species,
                            .packages = c("dbarts"),
                            .combine = rbind,
                            .errorhandling = "pass") %dopar% {
                              
                              y <- cover_matrix[[sp]]
                              
                              set.seed(which(modelable_species == sp))
                              folds <- sample(rep(1:5, length.out = length(y)))
                              
                              fold_metrics <- sapply(1:5, function(k) {
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
                                
                                obs <- y[test_idx]
                                rmse <- sqrt(mean((obs - pred)^2))
                                ss_res <- sum((obs - pred)^2)
                                ss_tot <- sum((obs - mean(obs))^2)
                                r2 <- 1 - ss_res / ss_tot
                                
                                c(rmse = rmse, r2 = r2)
                              })
                              
                              data.frame(species = sp, 
                                         cv_rmse = mean(fold_metrics["rmse", ]), 
                                         cv_r2 = mean(fold_metrics["r2", ]))
                            }

stopCluster(cl)

cv_results_cover <- cv_results_cover |> arrange(desc(cv_r2))
cv_results_cover |> tibble::as_tibble() |> print(n = Inf)

#### 0 - inflation #############################################################

zero_inflation <- cover_matrix |>
  dplyr::select(all_of(modelable_species)) |>
  summarise(across(everything(), ~ mean(.x == 0))) |>
  pivot_longer(everything(), names_to = "species", values_to = "prop_zero") |>
  arrange(desc(prop_zero))

zero_inflation |> print(n = Inf)
