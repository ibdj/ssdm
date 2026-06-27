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


#Step 2 — run cLHS for 50 points, then enforce ≥10 m spacing among them.

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

#Step 3 — enforce ≥10 m among the 50 selected points, dropping the lower-coverage-value one from any too-close pair.
# coverage-contribution proxy: cLHS objective ranks samples by selection order
sel$clhs_order <- seq_len(nrow(sel))   # earlier = pulled first by annealing

sel_sf <- st_as_sf(sel, coords = c("x", "y"), crs = 32622, remove = FALSE)

# pairwise distances among selected points
d <- st_distance(sel_sf) |> units::drop_units()
diag(d) <- Inf

# greedily drop the later-ordered point from any pair < 10 m
drop <- integer(0)
for (i in order(sel$clhs_order)) {
  if (i %in% drop) next
  too_close <- which(d[i, ] < 10)
  too_close <- setdiff(too_close, drop)
  too_close <- too_close[sel$clhs_order[too_close] > sel$clhs_order[i]]
  drop <- c(drop, too_close)
}

sel_keep <- if (length(drop) > 0) sel[-drop, ] else sel
nrow(sel_keep)   # should now be 50
length(drop)     # how many removed

#Step 4 — attach a coverage score and export for QGIS + handheld.Step 4 — attach a coverage score and export for QGIS + handheld.

# coverage score: lower clhs_order = pulled earlier = higher priority
sel_keep$priority <- rank(sel_keep$clhs_order)   # 1 = highest
sel_keep$plot_id  <- sprintf("NEW%03d", sel_keep$priority)

out_sf <- st_as_sf(sel_keep, coords = c("x", "y"), crs = 32622)

# GeoPackage for QGIS (keeps UTM + all attributes)
st_write(out_sf, "new_plots.gpkg", delete_dsn = TRUE)

# CSV in lat/lon for handheld GPS
ll <- st_transform(out_sf, 4326)
data.frame(
  plot_id  = ll$plot_id,
  priority = ll$priority,
  lat      = st_coordinates(ll)[, 2],
  lon      = st_coordinates(ll)[, 1],
  elevation = round(ll$elevation, 1),
  hli      = round(ll$hli, 2),
  ndvi     = round(ll$ndvi, 2)
) |> write.csv("new_plots.csv", row.names = FALSE)


e <- ext(rbind(vect(plots_sf), vect(out_sf))) * 1.1
plot(rast_ndvi_proc, ext = e)
plot(plots_sf, add = TRUE, col = "black", pch = 1)        # existing
plot(out_sf, add = TRUE, col = "red", pch = 19) 

