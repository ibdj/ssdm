#### packages ##################################################################

library(tidyverse)
library(vegan)
library(janitor)
library(terra)
library(sf)

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
