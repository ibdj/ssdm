#### packages ##################################################################

library(tidyverse)
library(vegan)
library(janitor)
library(terra)
library(sf)
library(car)
library(caret)

#### functions #################################################################

process_rast <- function(r, ref = ref_rast) {
  aoi_src <- project(aoi_masked, crs(r))   # AOI in the raster's own CRS
  r |>
    crop(aoi_src) |>                        # shrink BEFORE projecting
    project("EPSG:32622") |>
    mask(aoi_masked) |>
    resample(ref)
}

#### plots and aoi #############################################################

plots_sf <- readRDS("~/OneDrive - Aarhus universitet/MappingPlants/02 Modelling future changes/ssdm/data/plots_sf.rds")

aoi_plots <- plots_sf |>
  st_bbox() |>
  st_as_sfc() |>
  st_buffer(50) |>
  vect()  # convert to terra format for cropping

aoi_raw <- buffer(aoi_plots, width = 500)

plot(aoi_plots)

gpkg <- "data/aoi_masked.gpkg"

# see what layers the gpkg contains
vector_layers(gpkg)
aoi_masked <- vect(gpkg) 

#removing corner part
parts <- disagg(aoi_masked)
aoi_masked <- parts[which.max(expanse(parts)), ]

crs(aoi_masked, describe = TRUE)$code

#added 10 m buffer to the aoi
aoi_masked <- buffer(aoi_masked, width = 10)   # 30 m outward; units = metres (UTM)

plot(aoi_masked) #add = TRUE
plot(aoi_raw)

#### raster import #############################################################

rast_dem        <- rast("data/elevation_arcticdem-30_32622.tif") |> crop(aoi_masked)
rast_ndvi       <- rast("data/ndvi_export_2025.tif") |> crop(aoi_masked)
rast_ndwi       <- rast("data/ndwi.tif") 
rast_snowfree   <- rast("data/snow_free_days.tif")

rast_slope      <- terrain(rast_dem, v = "slope", unit = "degrees") |> crop(aoi_masked)

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

rast_temp_proc       <- temp_rast |> process_rast()

sapply(list(rast_dem_proc, 
            rast_ndvi_proc, 
            rast_ndwi_proc, 
            rast_snowfree_proc,
            rast_slope_proc,
            rast_hli_proc
            #rast_temp_proc
),
function(r) crs(r, describe = TRUE)$code)

plot(trim(rast_snowfree_proc))

#### raster solar radiation / heat load index ##################################

# hli needs to be calculated on the uncropped dem because there will be na cells in hli otherwise

hli <- spatialEco::hli(rast_dem)   # accepts a terra SpatRaster in recent versions
names(hli) <- "hli"
plot(hli)
plot(aoi_masked, add = TRUE) #add = TRUE

rast_hli_proc        <- hli |> process_rast()

#### raster temp (interpolation) ###############################################

# --- Temperature interpolation: predictor selection ---
# Goal: pick a small, physically sensible set of spatial predictors to
# interpolate mean growing-season soil temperature across the AOI.

# --- import the tms data and metrics ---
tms_all_sf <- readRDS("~/OneDrive - Aarhus universitet/MappingPlants/02 Modelling future changes/ssdm/data/tms_all_sf.rds") |> 
  mutate(serial_number = as.character(serial_number))
tms_metrics <- readRDS("~/OneDrive - Aarhus universitet/MappingPlants/02 Modelling future changes/ssdm/data/tms_metrics.rds")

names(tms_all_sf)

common <- intersect(names(tms_all_sf), names(tms_metrics))
print(common)   # expect serial (+ plot?), no value columns
tms_sf <- dplyr::left_join(tms_all_sf, tms_metrics, by = common)

# definting prediction stack 
pred_stack <- c(rast_dem_proc, rast_slope_proc, rast_ndwi_proc,
                rast_snowfree_proc, rast_hli_proc)
names(pred_stack) <- c("elevation", "slope", "ndwi", "snowfree", "hli")

tms_sf <- tms_sf |>
  dplyr::bind_cols(terra::extract(pred_stack, terra::vect(tms_sf), ID = FALSE))

summary(tms_sf)   # check: no NAs in the five new columns


#checking for multicolenearity
cand <- c("elevation", "hli", "ndwi", "snowfree", "slope")

# Correlation matrix among candidates: flags redundant predictors
# (e.g. ndvi/ndwi). Informs interpretation, not inclusion by itself.
round(cor(sf::st_drop_geometry(tms_sf)[, cand]), 2)

# VIF on the full candidate model: quantifies multicollinearity.
# Rule of thumb: >5 worth a look, >10 severe. Low VIF = safe to include.
full <- lm(summer_t ~ elevation + hli + ndwi + snowfree + slope,
           data = tms_sf)
vif(full)

# --- Best-subset selection by cross-validated prediction error ---
# Fit every non-empty subset of candidates; rank by LOOCV-RMSE.
# Out-of-sample error (not in-sample R2) is the selection criterion.

ctrl <- trainControl(method = "LOOCV", allowParallel = FALSE)

preds <- c("elevation", "hli", "ndwi", "snowfree", "slope")

# all non-empty subsets of the candidate predictors
combos <- unlist(
  lapply(seq_along(preds), \(k) combn(preds, k, simplify = FALSE)),
  recursive = FALSE
)

# LOOCV-RMSE for each candidate model
results <- sapply(combos, function(vars) {
  f <- reformulate(vars, response = "summer_t")
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
temp_lm <- lm(summer_t ~ hli + ndwi + slope, data = tms_sf)
summary(temp_lm)
par(mfrow = c(2, 2)); plot(temp_lm); par(mfrow = c(1, 1))

# interpolation with the model
rast_temp <- terra::predict(pred_stack, temp_lm)
names(rast_temp) <- "temp"
plot(trim(rast_temp))

summary(values(rast_temp), na.rm = TRUE)
plot(trim(rast_temp > 12.7))   # where the loggers' observed max is exceeded

