# ============================================================
# HMSC - Joint Species Distribution Model
# ============================================================

# Packages
library(Hmsc)
library(tidyverse)
library(terra)
library(sf)
library(corrplot)
library(parallel)

# ============================================================
# Load data from previous script outputs
# ============================================================

# Plot data
abiotic_plot    <- readRDS("data/abiotic_plot.rds")
species_matrix  <- readRDS("data/species_matrix.rds")
species_long    <- readRDS("data/species_long.rds")

# Rasters
dem_rast        <- rast("data/dem_crop.tif")
twi_rast        <- rast("data/twi_calculated.tif")
ndvi_rast       <- rast("data/ndvi_crop.tif")
slope_rast      <- rast("data/slope_crop.tif")
aspect_rast     <- rast("data/aspect_crop.tif")
temp_rast       <- rast("data/temp_predicted_rast.tif")

# ============================================================
# Prepare HMSC input data
# ============================================================

# Y matrix - presence/absence for modelable species
modelable_species <- c(
  "Coptis trifolia", "Huperzia selago", "Loiseleuria procumbens",
  "Salix arctophila", "Calamagrostis langsdorfii", "Juncus trifidus",
  "Ledum groenlandicum", "Polygonum viviparum", "Betula nana",
  "Lycopodium annotinum", "Salix herbacea", "Vaccinium uliginosum",
  "Salix glauca", "Deschampsia flexuosa", "Carex bigelowii",
  "Empetrum nigrum"
)

Y <- species_matrix |>
  dplyr::select(all_of(modelable_species)) |>
  mutate(across(everything(), ~ as.integer(. > 0))) |>
  as.matrix()

rownames(Y) <- species_matrix$plot_name

# X matrix - environmental predictors
X <- abiotic_plot |>
  dplyr::select(elevation, slope, aspect_sin, aspect_cos,
                twi, ndvi, temp_predicted) |>
  as.data.frame()

rownames(X) <- abiotic_plot$plot_name

# Spatial coordinates for random effect
coords <- abiotic_plot |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  st_transform(32622) |>
  st_coordinates() |>
  as.data.frame()

rownames(coords) <- abiotic_plot$plot_name

# Check dimensions align
dim(Y)
dim(X)
dim(coords)


# Check NAs in Y
colSums(is.na(Y))
rowSums(is.na(Y))

# Check NAs in X
colSums(is.na(X))
rowSums(is.na(X))

# Check NAs in coords
colSums(is.na(coords))

# ============================================================
# Set up HMSC model
# ============================================================

# Spatial random effect
studyDesign <- data.frame(plot = as.factor(rownames(Y)))
rL.spatial <- HmscRandomLevel(sData = coords)

# Formula for fixed effects
XFormula <- ~ elevation + slope + aspect_sin + aspect_cos + 
  twi + ndvi + temp_predicted

# Define model
m <- Hmsc(
  Y = Y,
  XData = X,
  XFormula = XFormula,
  studyDesign = studyDesign,
  ranLevels = list(plot = rL.spatial),
  distr = "probit"
)

# ============================================================
# MCMC sampling
# ============================================================

# Start with thin chain for testing
nChains   <- 4
thin      <- 10
samples   <- 1000
transient <- 500

set.seed(42)

m <- sampleMcmc(m,
                thin      = thin,
                samples   = samples,
                transient = transient,
                nChains   = nChains,
                verbose   = 500
)

# Save model immediately after sampling
saveRDS(m, "data/hmsc_model_pa.rds")

summary(m)
