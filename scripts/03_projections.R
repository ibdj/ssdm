#### packages ##################################################################

library(tidyverse)
library(terra)
library(dbarts) 
library(terra) 
library(raster)

#### packages ##################################################################

ssp_dirs <- c(ssp126 = "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/CMIP6/CMIP6 - Mean temperature (T) Change deg C - Long Term (2081-2100) SSP1-2.6 (rel. to 1995-2014) - June to August (32 models)",
              ssp245 = "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/CMIP6/CMIP6 - Mean temperature (T) Change deg C - Long Term (2081-2100) SSP2-4.5 (rel. to 1995-2014) - June to August (34 models)",
              ssp370 = "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/CMIP6/CMIP6 - Mean temperature (T) Change deg C - Long Term (2081-2100) SSP3-7.0 (rel. to 1995-2014) - June to August (30 models)",
              ssp585 = "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/CMIP6/CMIP6 - Mean temperature (T) Change deg C - Long Term (2081-2100) SSP5-8.5 (rel. to 1995-2014) - June to August (34 models)")

deltas <- sapply(ssp_dirs, \(d) {
  terra::extract(rast(file.path(d, "map.tif")),
                 matrix(c(-51.38079898859376, 64.13399461827079), ncol = 2))[1, 1]
})
deltas

#### setup #####################################################################
set.seed(42)

pred_stack <- rast("data/pred_stack.tif")
pred_names <- names(pred_stack)

modelable_species <- readRDS("data/species_frequency.rds") |>
  dplyr::filter(n_plots >= 10) |>
  dplyr::pull(taxon)

deltas <- c(ssp126 = 0.9428791, ssp245 = 1.9663557,
            ssp370 = 2.8462466, ssp585 = 3.8888474)

for (sc in names(deltas)) {
  
  stack_sc <- pred_stack
  stack_sc[["summer_t_int"]] <- stack_sc[["summer_t_int"]] + deltas[[sc]]
  
  stack_sc_r <- raster::stack(stack_sc)
  pred_df    <- as.data.frame(stack_sc_r, na.rm = FALSE)
  idx        <- complete.cases(pred_df)
  
  for (sp in modelable_species) {
    model <- readRDS(paste0("data/bart_pa_", gsub(" ", "_", sp), ".rds"))
    p     <- colMeans(dbarts:::predict.bart(model, newdata = pred_df[idx, ]))
    
    full  <- rep(NA_real_, nrow(pred_df)); full[idx] <- p
    r     <- raster::raster(stack_sc_r[[1]]); raster::values(r) <- full
    raster::writeRaster(r,
                        paste0("data/sdm_", gsub(" ", "_", sp), "_", sc, ".tif"), overwrite = TRUE)
  }
  cat(sc, "done\n")
}

#### change maps and richness stacks ###########################################
scenarios <- names(deltas)

# stacked richness (sum of PA probabilities) per time slice
richness <- list()

richness$present <- rast(lapply(modelable_species, \(sp)
                                rast(paste0("data/sdm_", gsub(" ", "_", sp), ".tif")))) |> sum()

for (sc in scenarios) {
  richness[[sc]] <- rast(lapply(modelable_species, \(sp)
                                rast(paste0("data/sdm_", gsub(" ", "_", sp), "_", sc, ".tif")))) |> sum()
  writeRaster(richness[[sc]] - richness$present,
              paste0("data/richness_change_", sc, ".tif"), overwrite = TRUE)
}
writeRaster(richness$present, "data/richness_present.tif", overwrite = TRUE)

# per-species change, headline scenario example
plot(trim(rast(paste0("data/sdm_Betula_nana_ssp585.tif")) -
            rast("data/sdm_Betula_nana.tif")))
