# 06_richness_future.R

library(terra)
library(tidyverse)

modelable_species <- readRDS("data/modelable_species.rds")
scenarios <- c("current", "ssp126", "ssp245", "ssp370", "ssp585")
deltas <- c(current = 0, ssp126 = 0.9428791, ssp245 = 1.9663557,
            ssp370 = 2.8462466, ssp585 = 3.8888474)

out_files <- c(rast = "data/richness_future.tif",
               post = "data/richness_posterior.rds",
               pix  = "data/richness_change_pixels.rds")

#### species change vs borealness ##############################################

species_change <- readRDS("data/boreal_sp_change_future.rds")
stopifnot(sum(is.na(species_change$B_std)) == 0)

ggplot(species_change, aes(B_std, delta_p)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40") +
  geom_point() +
  ggrepel::geom_text_repel(
    data = dplyr::filter(species_change, abs(delta_p) > 0.05),
    aes(label = taxon), size = 2.5, fontface = "italic"
  ) +
  facet_wrap(~scenario) +
  labs(x = "Borealness (B_std)",
       y = "Change in mean occurrence probability (NDVI-weighted)") +
  theme_minimal()

#### richness maps (from posterior-mean SDM rasters) ###########################

if (!file.exists(out_files["rast"])) {
  richness <- purrr::map(scenarios, \(s) {
    f <- file.path("data", paste0("sdm_", gsub(" ", "_", modelable_species),
                                  if (s == "current") "" else paste0("_", s), ".tif"))
    sum(terra::rast(f))
  }) |> terra::rast()
  names(richness) <- scenarios
  terra::writeRaster(richness, out_files["rast"], overwrite = TRUE)
}

richness    <- terra::rast(out_files["rast"])
rich_change <- richness[[-1]] - richness[["current"]]

terra::global(richness[["current"]], c("min", "max", "mean"), na.rm = TRUE)

#### posterior richness change (only if outputs are missing) ###################

if (!all(file.exists(out_files[c("post", "pix")]))) {
  
  pred_stack <- terra::rast("data/pred_stack.tif")
  pred_df    <- terra::as.data.frame(pred_stack, na.rm = FALSE)
  idx        <- complete.cases(pred_df)
  
  draws <- seq(4, 2000, by = 4)   # 500 thinned draws, as for borealness
  pos   <- split(seq_len(sum(idx)), ceiling(seq_len(sum(idx)) / 10000))
  
  models <- purrr::map(purrr::set_names(modelable_species), \(sp)
                       readRDS(paste0("data/bart_pa_", gsub(" ", "_", sp), ".rds")))
  
  # saved models must keep their tree state (not predict 0.5 everywhere)
  stopifnot(purrr::map_lgl(models, \(m)
                           diff(range(dbarts:::predict.bart(m, newdata = pred_df[idx, ][1:5, ]))) > 0))
  
  rich_post <- function(delta) {
    df <- pred_df[idx, ]
    df$summer_t_int <- df$summer_t_int + delta
    out <- matrix(0, length(draws), nrow(df))
    for (m in models) {
      for (ps in pos) {
        out[, ps] <- out[, ps] +
          dbarts:::predict.bart(m, newdata = df[ps, ])[draws, , drop = FALSE]
      }
    }
    out
  }
  
  rich_cur <- rich_post(deltas["current"])
  
  res <- purrr::map(purrr::set_names(names(deltas)[-1]), \(s) {
    d <- rich_post(deltas[s]) - rich_cur
    list(post = tibble::tibble(scenario = s,
                               draw     = seq_along(draws),
                               cur      = rowMeans(rich_cur),
                               change   = rowMeans(d),
                               share_up = rowMeans(d > 0)),
         pix  = matrixStats::colMedians(d))
  })
  
  saveRDS(purrr::map_dfr(res, "post"), out_files["post"])
  
  pix_med <- purrr::map(res, "pix")
  saveRDS(tibble::tibble(scenario = rep(names(pix_med), each = sum(idx)),
                         change   = unlist(pix_med, use.names = FALSE)),
          out_files["pix"])
}

#### load results ##############################################################

rich_res    <- readRDS(out_files["post"])
rich_pix_df <- readRDS(out_files["pix"])

#### credible intervals ########################################################

rich_ci <- rich_res |>
  dplyr::summarise(med    = median(change),
                   lo     = quantile(change, 0.025),
                   hi     = quantile(change, 0.975),
                   p_gain = mean(change > 0),
                   .by = scenario)
rich_ci

median(rich_res$cur)   # current mean richness across the study area

rich_pix_df |> dplyr::summarise(share_up = mean(change > 0), .by = scenario)

#### histogram of per-pixel change #############################################

ggplot(rich_pix_df, aes(x = change, fill = scenario)) +
  geom_histogram(bins = 60, position = "identity", alpha = 0.4) +
  geom_vline(xintercept = 0, colour = "grey30") +
  geom_vline(data = rich_ci, aes(xintercept = med, colour = scenario), linewidth = 0.8) +
  geom_vline(data = rich_ci, aes(xintercept = lo, colour = scenario), linetype = "dashed") +
  geom_vline(data = rich_ci, aes(xintercept = hi, colour = scenario), linetype = "dashed") +
  labs(x = "Change in stacked species richness (posterior median)",
       y = "Number of pixels",
       fill = NULL, colour = NULL) +
  theme_minimal()
