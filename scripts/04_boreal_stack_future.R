# 04_boreal_stack_future.R

library(tidyverse)
library(terra)

# matching names


modelable_species <- readRDS("data/modelable_species.rds")
idx_modelable     <- readRDS("data/idx_modelable.rds")

scenarios <- c("current", "ssp126", "ssp245", "ssp370", "ssp585")

files <- tidyr::expand_grid(taxon = modelable_species, scenario = scenarios) |>
  dplyr::mutate(file = file.path("data", paste0("sdm_", gsub(" ", "_", taxon),
                                                dplyr::if_else(scenario == "current", "", paste0("_", scenario)),
                                                ".tif")))

sum(!file.exists(files$file))   # 0 expected
files$file[!file.exists(files$file)]

# Step 2: computing the probability-weighted mean borealness (CWM) per pixel for each scenario, and save it

b <- idx_modelable$B_std[match(modelable_species, idx_modelable$taxon)]
stopifnot(!anyNA(b))

cwm <- purrr::map(scenarios, \(s) {
  r <- terra::rast(files$file[files$scenario == s])
  terra::app(r, \(p) sum(p * b) / sum(p))
}) |>
  terra::rast()
names(cwm) <- scenarios

terra::writeRaster(cwm, "data/boreal_cwm_future.tif", overwrite = TRUE)

plot(trim(cwm))
global(cwm, c("mean", "sd"), na.rm = TRUE)

# Step 3: map and summarise change from current

cwm_change <- cwm[[-1]] - cwm[["current"]]
names(cwm_change) <- scenarios[-1]

terra::writeRaster(cwm_change, "data/boreal_cwm_change_future.tif", overwrite = TRUE)

plot(trim(cwm_change), col = hcl.colors(50, "Blue-Red 2"),
     range = c(-1, 1) * max(abs(global(cwm_change, "range", na.rm = TRUE))))

global(cwm_change, \(x) quantile(x, c(0.05, 0.5, 0.95), na.rm = TRUE))
global(cwm_change > 0, "mean", na.rm = TRUE)   # share of pixels becoming more boreal
