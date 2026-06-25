library(clhs)

# --- Step 1: build the candidate pool ---

# stack the three placement layers (already masked to aoi_masked, EPSG:32622)
place_stack <- c(rast_dem_proc, rast_hli_proc, rast_ndvi_proc)
names(place_stack) <- c("elevation", "hli", "ndvi")

# every non-NA cell -> data frame with cell coordinates
cand <- as.data.frame(place_stack, xy = TRUE, na.rm = TRUE)
nrow(cand)   # eligible cells before exclusion

# drop candidate cells within 10 m of any existing plot
cand_sf <- st_as_sf(cand, coords = c("x", "y"), crs = 32622, remove = FALSE)

nearest_dist <- st_distance(cand_sf, plots_sf) |>
  apply(1, min)                              # min distance to an existing plot (metres)

cand <- cand[nearest_dist >= 10, ]           # keep cells >=10 m from existing plots
nrow(cand)   # eligible cells after exclusion

set.seed(42)   # reproducible; recorded in script

res <- clhs(
  cand[, c("elevation", "hli", "ndvi")],   # only the placement vars
  size = 50,
  iter = 10000,
  progress = FALSE,
  simple = FALSE
)

# pull selected rows (with coordinates) back out
sel <- cand[res$index_samples, ]
nrow(sel)   # should be 50