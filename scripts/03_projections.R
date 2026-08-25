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