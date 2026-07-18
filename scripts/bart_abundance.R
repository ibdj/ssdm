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

#### cover matrix ##############################################################
# Build cover matrix (replaces pa_matrix)
cover_matrix <- species_matrix |>
  dplyr::select(plot_name, all_of(modelable_species)) |>
  left_join(
    mp_abiotic |>
      dplyr::select(plot_name, elevation, slope, hli, ndwi, temp, snowfree),
    by = "plot_name"
  )

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


