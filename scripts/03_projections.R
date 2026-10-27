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

# Plot change for first scenario
plot(trim(rast(paste0("data/richness_present.tif"))))
plot(trim(rast(paste0("data/richness_change_", scenarios[1], ".tif"))))
plot(trim(rast(paste0("data/richness_change_", scenarios[2], ".tif"))))
plot(trim(rast(paste0("data/richness_change_", scenarios[3], ".tif"))))


# per-species change, headline scenario example
plot(trim(rast(paste0("data/sdm_Betula_nana_ssp585.tif")) -
            rast("data/sdm_Betula_nana.tif")))

plot(trim(rast(paste0("data/sdm_Betula_nana_ssp126.tif")) -
            rast("data/sdm_Betula_nana.tif")))

sapply(names(deltas), \(sc)
       global(richness[[sc]] - richness$present, "mean", na.rm = TRUE)[1, 1])


#### pr species change #########

sp_change <- map_dfr(modelable_species, function(sp) {
  pres <- rast(paste0("data/sdm_", gsub(" ", "_", sp), ".tif"))
  map_dfr(names(deltas), function(sc) {
    fut <- rast(paste0("data/sdm_", gsub(" ", "_", sp), "_", sc, ".tif"))
    tibble(species = sp, scenario = sc,
           mean_change = global(fut - pres, "mean", na.rm = TRUE)[1, 1])
  })
})

sp_change |>
  pivot_wider(names_from = scenario, values_from = mean_change) |>
  arrange(desc(ssp585)) |>
  print(n = 17)

#### partial dependence: summer_t_int ##########################################
x_train <- readRDS("data/rs.rds") |> sf::st_drop_geometry() |>
  dplyr::select(dplyr::all_of(pred_names)) |> as.data.frame()

t_grid <- seq(min(x_train$summer_t_int), max(x_train$summer_t_int) + 3.9,
              length.out = 40)   # extend into projection range

pd <- map_dfr(modelable_species, function(sp) {
  model <- readRDS(paste0("data/bart_pa_", gsub(" ", "_", sp), ".rds"))
  map_dfr(t_grid, function(tv) {
    x_pd <- x_train
    x_pd$summer_t_int <- tv
    tibble(species = sp, summer_t_int = tv,
           p = mean(colMeans(dbarts:::predict.bart(model, newdata = x_pd))))
  })
})

ggplot(pd, aes(summer_t_int, p)) +
  geom_line() +
  geom_vline(xintercept = max(x_train$summer_t_int), linetype = 2, colour = "grey50") +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = "Interpolated summer soil temperature (°C)",
       y = "Mean predicted P(presence)") +
  theme_minimal() +
  theme(strip.text = element_text(face = "italic"))

cover_matrix |>   # or pa_matrix
  tidyr::pivot_longer(dplyr::all_of(modelable_species), names_to = "species", values_to = "pa") |>
  dplyr::filter(pa > 0) |>
  dplyr::summarise(ndwi_mean = mean(ndwi), .by = species) |> arrange(desc(ndwi_mean))
  dplyr::arrange(dplyr::desc(ndwi_mean))  

  #dplyr::summarise(t_mean = mean(summer_t_int), n = dplyr::n(), .by = species) |>
  #dplyr::arrange(dplyr::desc(t_mean))
# richness
  
  # Mean richness change per scenario
  richness_summary <- map_dfr(scenarios, function(sc) {
    change <- rast(paste0("data/richness_change_", sc, ".tif"))
    tibble(
      scenario = sc,
      mean_change = mean(values(change), na.rm = TRUE),
      sd_change   = sd(values(change), na.rm = TRUE),
      pct_increase = mean(values(change) > 0, na.rm = TRUE) * 100,
      pct_decrease = mean(values(change) < 0, na.rm = TRUE) * 100
    )
  })
  
  print(richness_summary)

  richness_summary |>
    pivot_longer(cols = c(pct_increase, pct_decrease),
                 names_to = "direction",
                 values_to = "pct") |>
    mutate(direction = factor(direction, 
                              levels = c("pct_increase", "pct_decrease"),
                              labels = c("Increase", "Decrease"))) |>
    ggplot(aes(x = scenario, y = pct, fill = direction)) +
    geom_bar(stat = "identity") +
    geom_text(aes(x = scenario, label = round(mean_change, 3)), 
              data = richness_summary,
              y = 50, colour = "white", fontface = "bold",
              inherit.aes = FALSE)
    scale_fill_manual(values = c("Increase" = "#2E75B6", "Decrease" = "#C00000")) +
    labs(x = "Scenario", y = "% of pixels", fill = NULL,
         title = "Richness change direction per scenario",
         subtitle = "Numbers show mean richness change") +
    theme_minimal()
    
    # numero of clusters
    
    species_matrix <- readRDS("data/species_matrix.rds")
    
    # Hellinger transform for clustering (handles double zeros well)
    sp_hell <- decostand(species_matrix |> dplyr::select(-plot_name), 
                         method = "hellinger")
    
    # Hierarchical clustering
    dist_mat <- vegdist(sp_hell, method = "bray")
    clust <- hclust(dist_mat, method = "ward.D2")
    plot(clust)
    
    # Using cascadeKM from vegan
    casc <- cascadeKM(sp_hell, inf.gr = 2, sup.gr = 10, iter = 100)
    plot(casc)
    casc$results  # shows SSE and Calinski criterion per k
    
    # environment in clusters
    
    clusters <- cutree(clust, k = 5)
    saveRDS(clusters, "data/clusters_k5.rds")
    
    plot_clusters <- tibble(
      plot_name = species_matrix$plot_name,
      cluster = factor(clusters)
    )
    # Check environmental separation
    plot_clusters |>
      left_join(readRDS("data/abiotic_plot.rds"), by = "plot_name") |>
      group_by(cluster) |>
      summarise(
        n = n(),
        mean_temp = mean(temp_mean_tms, na.rm = TRUE),
        mean_elev = mean(elevation, na.rm = TRUE),
        mean_snowfree = mean(snowfree, na.rm = TRUE)
      )
  # bart on clusters
    
    rs <- readRDS("data/rs.rds") |>
      st_drop_geometry() |>
      rename(plot_name = plot)
    
    pa_clusters <- plot_clusters |>
      left_join(rs, by = "plot_name")
    
    # Verify predictors are all there
    pa_clusters |> 
      dplyr::select(all_of(names(pred_stack))) |>
      summarise(across(everything(), ~sum(is.na(.x))))
    
    # Add clusters to abiotic data
    pa_clusters <- plot_clusters |>
      left_join(readRDS("data/abiotic_plot.rds"), by = "plot_name")
    
    # Load pred stack
    pred_stack <- rast("data/pred_stack.tif")
    pred_stack_r <- raster::stack(pred_stack)
    pred_df <- as.data.frame(pred_stack_r, na.rm = FALSE)
    idx <- complete.cases(pred_df)
    
    # Fit and project one model per cluster
    cluster_levels <- levels(pa_clusters$cluster)
    
    for (cl in cluster_levels) {
      cat("Fitting cluster:", cl, "\n")
      
      y <- as.integer(pa_clusters$cluster == cl)
      
      train_data <- cbind(
        pa_clusters |> dplyr::select(all_of(names(pred_stack))),
        y = y
      ) |> as.data.frame()
      
      model <- dbarts::bart2(y ~ ., data = train_data, 
                             keepTrees = TRUE, seed = as.integer(cl))
      
      # Present projection
      p <- colMeans(dbarts:::predict.bart(model, newdata = pred_df[idx, ]))
      full <- rep(NA_real_, nrow(pred_df)); full[idx] <- p
      r <- raster::raster(pred_stack_r[[1]]); raster::values(r) <- full
      raster::writeRaster(r, paste0("data/community_", cl, "_present.tif"), 
                          overwrite = TRUE)
      
      # Scenario projections
      deltas <- c(ssp126 = 0.9428791, ssp245 = 1.9663557,
                  ssp370 = 2.8462466, ssp585 = 3.8888474)
      
      for (sc in names(deltas)) {
        stack_sc <- pred_stack
        stack_sc[["summer_t_int"]] <- stack_sc[["summer_t_int"]] + deltas[[sc]]
        stack_sc_r <- raster::stack(stack_sc)
        pred_df_sc <- as.data.frame(stack_sc_r, na.rm = FALSE)
        
        p_sc <- colMeans(dbarts:::predict.bart(model, newdata = pred_df_sc[idx, ]))
        full_sc <- rep(NA_real_, nrow(pred_df_sc)); full_sc[idx] <- p_sc
        r_sc <- raster::raster(pred_stack_r[[1]]); raster::values(r_sc) <- full_sc
        raster::writeRaster(r_sc, paste0("data/community_", cl, "_", sc, ".tif"),
                            overwrite = TRUE)
      }
    }
    
    par(mfrow = c(2, 3))
    for (cl in 1:5) {
      plot(trim(rast(paste0("data/community_", cl, "_present.tif"))),
           main = paste("Community", cl))
    }
    par(mfrow = c(1, 1))
    
    # Species summary per cluster
    species_summary <- species_matrix |>
      dplyr::select(-plot_name) |>
      mutate(cluster = clusters) |>
      group_by(cluster) |>
      summarise(across(everything(), mean)) |>
      pivot_longer(-cluster, names_to = "species", values_to = "mean_cover") |>
      group_by(cluster) |>
      slice_max(mean_cover, n = 3) |>
      summarise(dominant_species = paste(species, collapse = ", "))
    
    # Abiotic summary per cluster
    abiotic_summary <- pa_clusters |>
      group_by(cluster) |>
      summarise(
        n = n(),
        temp = round(mean(summer_t_int, na.rm = TRUE), 1),
        elevation = round(mean(elevation, na.rm = TRUE), 0),
        snowfree = round(mean(snowfree, na.rm = TRUE), 0),
        hli = round(mean(hli, na.rm = TRUE), 2)
      )
    
    # Combined table
    community_table <- abiotic_summary |>
      left_join(species_summary, by = "cluster") |>
      dplyr::select(cluster, n, dominant_species, temp, elevation, snowfree, hli)
    
    community_table
    
    # mean B of the "community" now vs projected
    B_now    <- weighted.mean(idx_modelable$B_std, w = df$area_now)
    B_future <- weighted.mean(idx_modelable$B_std, w = df$area_future)
    B_future - B_now     # > 0 = the area as a whole borealizes
    