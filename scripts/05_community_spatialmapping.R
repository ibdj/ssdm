library(terra)

# Load data
species_matrix <- readRDS("data/species_matrix_pa.rds")
clusters       <- readRDS("data/clusters_k5.rds")
modelable      <- readRDS("data/modelable_species.rds")

# Cluster centroids (mean PA per species per cluster)
sp_vals <- species_matrix[, colnames(species_matrix) != "plot_name"]

centroids <- do.call(rbind, lapply(1:5, function(k) {
  colMeans(sp_vals[clusters == k, , drop = FALSE])
}))
rownames(centroids) <- paste0("C", 1:5)

# Stack present rasters (only modelable species)
raster_files <- file.path("data", paste0("sdm_", gsub(" ", "_", modelable), ".tif"))
r_stack <- rast(raster_files)
names(r_stack) <- modelable

# Convert stack to matrix of pixel values
px <- as.matrix(r_stack)

# Euclidean distance from each pixel to each centroid
# centroids cols must match raster layer order
centroid_sub <- centroids[, modelable, drop = FALSE]

dist_to_centroid <- apply(px, 1, function(p) {
  sqrt(rowSums((centroid_sub - p)^2))
})
# dist_to_centroid is 5 x npixels — transpose
dist_to_centroid <- t(dist_to_centroid)

# Assign nearest centroid
assigned <- apply(dist_to_centroid, 1, which.min)

# Min distance for each pixel (how far from nearest centroid)
min_dist <- apply(dist_to_centroid, 1, min)

# Transitional threshold: e.g. 75th percentile of min distances
threshold <- quantile(min_dist, 0.75, na.rm = TRUE)

# 0 = transitional, 1-5 = community
community_assigned <- ifelse(min_dist > threshold, 0L, as.integer(assigned))

# Write back to raster
r_community <- r_stack[[1]]
values(r_community) <- community_assigned
names(r_community) <- "community"

writeRaster(r_community, "data/community_present.tif", overwrite = TRUE)

r_cropped <- trim(r_community)
plot(r_cropped, col = c("grey80", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
     breaks = -0.5:5.5, legend = TRUE,
     main = "Community assignment (0 = transitional)")

table(values(r_community))


######## community #############

# Prevalence of each species within each cluster
prev_by_cluster <- do.call(rbind, lapply(1:5, function(k) {
  colMeans(sp_vals[clusters == k, , drop = FALSE])
}))
rownames(prev_by_cluster) <- paste0("C", 1:5)

# Top 4 species per community
indicators <- lapply(1:5, function(k) {
  row <- prev_by_cluster[k, ]
  top4 <- names(sort(row, decreasing = TRUE)[1:4])
  data.frame(species = top4, threshold = row[top4])
})
names(indicators) <- paste0("C", 1:5)

indicators

# For each community, check if all 4 indicators exceed threshold
px_df <- as.data.frame(px)

community_manual <- sapply(1:nrow(px_df), function(i) {
  pixel <- as.numeric(px_df[i, ])
  names(pixel) <- colnames(px_df)
  
  if (any(is.na(pixel))) return(NA_integer_)
  
  matches <- sapply(1:5, function(k) {
    ind <- indicators[[k]]
    ind <- ind[ind$species %in% names(pixel), ]
    all(pixel[ind$species] >= ind$threshold)
  })
  
  if (sum(matches) == 1) as.integer(which(matches))
  else 0L
})

r_community_manual <- r_stack[[1]]
values(r_community_manual) <- as.integer(community_manual)

plot(trim(r_community_manual), 
     col = c("grey80", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
     breaks = -0.5:5.5, legend = TRUE,
     main = "Community assignment - manual thresholds (0 = transitional)")

library(indicspecies)
species_matrix <- readRDS("data/species_matrix.rds")
sp_vals <- species_matrix[, colnames(species_matrix) != "plot_name"]

indval <- multipatt(as.data.frame(sp_vals), clusters, func = "IndVal.g")
summary(indval, indvalcomp = TRUE)

indicators_indval2 <- list(
  C1 = c("Vaccinium uliginosum", "Betula nana"),
  C2 = c("Salix glauca", "Calamagrostis langsdorfii"),
  C3 = NULL,
  C4 = c("Juncus trifidus"),
  C5 = NULL
)

# Thresholds: prevalence of each indicator in its own cluster
thresholds_indval2 <- lapply(thresholds_indval, function(df) {
  if (is.null(df)) return(NULL)
  df$threshold <- df$threshold / 100
  df
})
names(thresholds_indval2) <- paste0("C", 1:5)

community_manual2 <- sapply(1:nrow(px_df), function(i) {
  pixel <- as.numeric(px_df[i, ])
  names(pixel) <- colnames(px_df)
  
  if (any(is.na(pixel))) return(NA_integer_)
  
  matches <- sapply(1:5, function(k) {
    ind <- thresholds_indval2[[k]]
    if (is.null(ind)) return(FALSE)
    vals <- pixel[ind$species]
    if (any(is.na(vals))) return(FALSE)
    all(vals >= ind$threshold)
  })
  if (sum(matches) == 1) as.integer(which(matches))
  else 0L
})

r_community_manual2 <- r_stack[[1]]


values(r_community_manual2) <- as.integer(community_manual2)

plot(trim(r_community_manual2),
     col = c("grey80", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
     breaks = -0.5:5.5, legend = TRUE,
     main = "Community assignment - IndVal thresholds (0 = transitional)")

library(randomForest)

rs <- readRDS("data/rs.rds")

train_df <- data.frame(
  plot_name = species_matrix[, "plot_name"],
  sp_vals
) |>
  dplyr::left_join(
    rs |> sf::st_drop_geometry() |> dplyr::select(plot, elevation, hli, snowfree),
    by = c("plot_name" = "plot")
  ) |>
  dplyr::mutate(community = factor(clusters)) |>
  dplyr::select(-plot_name)

set.seed(42)
rf_model <- randomForest(
  community ~ .,
  data = train_df |> tidyr::drop_na(),
  ntree = 500,
  importance = TRUE
)

print(rf_model)

modelable_fixed <- gsub("Bistorta vivipara", "Polygonum viviparum", modelable)
model_cols <- make.names(modelable_fixed)

names(train_df) <- make.names(names(train_df))

train_df2 <- train_df |> 
  dplyr::select(dplyr::all_of(c(model_cols, "elevation", "hli", "snowfree", "community")))

set.seed(42)
rf_model2 <- randomForest(
  community ~ .,
  data = train_df2 |> tidyr::drop_na(),
  ntree = 500,
  importance = TRUE
)

print(rf_model2)

# Rename r_full layers to match model
names(r_full) <- make.names(c(gsub(" ", "_", modelable_fixed) |> gsub("_", ".", x=_), "elevation", "hli", "snowfree"))

r_rf <- predict(r_full, rf_model2, type = "response", na.rm = TRUE)

plot(trim(r_rf),
     col = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
     main = "Community RF classification")

train_df <- data.frame(
  plot_name = species_matrix[, "plot_name"],
  sp_vals
) |>
  dplyr::left_join(
    abiotic[, c("plot_name", "elevation", "snowfree", 
                "veg_height_ave", "ndvi", "bare_ground_bb")],
    by = "plot_name"
  ) |>
  dplyr::left_join(
    rs |> sf::st_drop_geometry() |> dplyr::select(plot, hli),
    by = c("plot_name" = "plot")
  ) |>
  dplyr::mutate(community = factor(clusters)) |>
  dplyr::select(-plot_name)

names(train_df) <- make.names(names(train_df))

train_df3 <- train_df |>
  dplyr::select(dplyr::all_of(c(model_cols, "elevation", "hli", "snowfree",
                                "veg_height_ave", "ndvi", "bare_ground_bb", "community")))

set.seed(42)
rf_model3 <- randomForest(
  community ~ .,
  data = train_df3 |> tidyr::drop_na(),
  ntree = 500,
  importance = TRUE
)

print(rf_model3)

train_df4 <- train_df |>
  dplyr::select(dplyr::all_of(c(model_cols, "elevation", "hli", "snowfree", "ndvi", "community")))

set.seed(42)
rf_model4 <- randomForest(
  community ~ .,
  data = train_df4 |> tidyr::drop_na(),
  ntree = 500,
  importance = TRUE
)

print(rf_model4)
