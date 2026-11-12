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

# load NDVI as a weight, check alignment, and compute NDVI-weighted means

ndvi <- terra::rast("data/rast_ndvi_proc.tif")
terra::compareGeom(ndvi, cwm)   # TRUE expected

w <- terra::clamp(ndvi, lower = 0, values = TRUE) |>
  terra::mask(cwm[[1]])

wmean <- \(x) terra::global(x * w, "sum", na.rm = TRUE)[[1]] /
  terra::global(w, "sum", na.rm = TRUE)[[1]]

tibble::tibble(
  scenario  = scenarios,
  cwm_unw   = terra::global(cwm, "mean", na.rm = TRUE)[[1]],
  cwm_ndvi  = purrr::map_dbl(scenarios, \(s) wmean(cwm[[s]]))
) |>
    dplyr::mutate(change_unw  = cwm_unw  - cwm_unw[1],
                  change_ndvi = cwm_ndvi - cwm_ndvi[1])

# species drivers of borealness changes
sp_change <- files |>
  dplyr::mutate(p_mean = purrr::map_dbl(file, \(f) wmean(terra::rast(f)))) |>
  dplyr::select(taxon, scenario, p_mean) |>
  tidyr::pivot_wider(names_from = scenario, values_from = p_mean) |>
  tidyr::pivot_longer(dplyr::starts_with("ssp"), names_to = "scenario", values_to = "p_future") |>
  dplyr::mutate(delta_p = p_future - current) |>
  dplyr::left_join(dplyr::select(idx_modelable, taxon, B_std), by = "taxon")

saveRDS(sp_change, "data/boreal_sp_change_future.rds")

ggplot(sp_change, aes(B_std, delta_p)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40") +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = taxon), size = 2.5) +
  facet_wrap(~scenario) +
  labs(x = "Borealness (B_std)", y = "Change in mean occurrence probability") +
  theme_minimal()

# 05_boreal_cwm_posterior.R
# Posterior credible intervals for AOI-level change in community borealness (CWM)

library(tidyverse)
library(terra)
library(dbarts)

#### setup #####################################################################

modelable_species <- readRDS("data/modelable_species.rds")
idx_modelable     <- readRDS("data/idx_modelable.rds")

b <- idx_modelable$B_std[match(modelable_species, idx_modelable$taxon)]
stopifnot(!anyNA(b))

pred_stack <- terra::rast("data/pred_stack.tif")
pred_names <- names(pred_stack)

deltas <- c(current = 0, ssp126 = 0.9428791, ssp245 = 1.9663557,
            ssp370 = 2.8462466, ssp585 = 3.8888474)

pred_df <- terra::as.data.frame(pred_stack, na.rm = FALSE)
idx     <- complete.cases(pred_df)
sum(idx)   # 63230 expected

#### refit BART models in session ##############################################
# saved bart_models.rds has lost its C++ tree state (predictions = 0.5)

pa_matrix <- readRDS("data/species_matrix_pa.rds") |>
  as.data.frame() |>
  tibble::rownames_to_column("plot") |>
  dplyr::select(plot, dplyr::all_of(modelable_species)) |>
  dplyr::left_join(sf::st_drop_geometry(readRDS("data/rs.rds")), by = "plot")

x_train <- pa_matrix |>
  dplyr::select(dplyr::all_of(pred_names)) |>
  as.data.frame()

bart_models <- purrr::map(purrr::set_names(modelable_species), \(sp)
                          dbarts::bart2(y ~ ., data = cbind(x_train, y = pa_matrix[[sp]]),
                                        keepTrees = TRUE,
                                        seed = which(modelable_species == sp),
                                        verbose = FALSE))

p_test <- dbarts:::predict.bart(bart_models[[1]], newdata = pred_df[idx, ][1:5, ])
range(p_test)   # must vary, not 0.5-0.5

#### posterior CWM per draw ####################################################

draws <- seq(4, 2000, by = 4)   # 500 thinned draws
pos   <- split(seq_len(sum(idx)), ceiling(seq_len(sum(idx)) / 10000))

ndvi <- terra::rast("data/rast_ndvi_proc.tif")
w <- pmax(terra::values(ndvi, mat = FALSE)[idx], 0)
w[is.na(w)] <- 0

cwm_post <- function(delta) {
  df <- pred_df[idx, ]
  df$summer_t_int <- df$summer_t_int + delta
  num <- den <- matrix(0, length(draws), nrow(df))
  for (i in seq_along(modelable_species)) {
    for (ps in pos) {
      p <- dbarts:::predict.bart(bart_models[[i]], newdata = df[ps, ])[draws, , drop = FALSE]
      num[, ps] <- num[, ps] + p * b[i]
      den[, ps] <- den[, ps] + p
    }
    cat(names(delta), modelable_species[i], "\n")
  }
  num / den
}

cwm_cur <- cwm_post(deltas["current"])

cwm_res <- purrr::map_dfr(names(deltas)[-1], \(s) {
  d <- cwm_post(deltas[s]) - cwm_cur
  tibble::tibble(scenario   = s,
                 draw       = seq_along(draws),
                 cur_unw    = rowMeans(cwm_cur),
                 delta_unw  = rowMeans(d),
                 delta_ndvi = as.vector(d %*% w) / sum(w),
                 share_up   = rowMeans(d > 0))
})

saveRDS(cwm_res, "data/boreal_cwm_posterior.rds")

#### summarise credible intervals ##############################################

sd(cwm_res$cur_unw[cwm_res$scenario == "ssp126"])

cwm_res |>
  dplyr::group_by(scenario) |>
  dplyr::summarise(
    dplyr::across(c(delta_unw, delta_ndvi, share_up),
                  list(med = median,
                       lo  = \(x) quantile(x, 0.025),
                       hi  = \(x) quantile(x, 0.975))),
    p_decline      = mean(delta_unw < 0),
    p_decline_ndvi = mean(delta_ndvi < 0)) |>
  print(width = Inf)

# histogram of pr pixcel change
pix_med <- purrr::map(purrr::set_names(names(deltas)[-1]), \(s) {
  d <- cwm_post(deltas[s]) - cwm_cur
  matrixStats::colMedians(d)
})

med_rast <- purrr::map(pix_med, \(v) {
  r <- terra::rast(pred_stack[[1]])
  r[which(idx)] <- v
  r
}) |>
  terra::rast()
names(med_rast) <- names(pix_med)
terra::writeRaster(med_rast, "data/boreal_cwm_change_median.tif", overwrite = TRUE)

pix_df <- tibble::tibble(scenario = rep(names(pix_med), each = sum(idx)),
                         change   = unlist(pix_med, use.names = FALSE))
saveRDS(pix_df, "data/boreal_cwm_change_pixels.rds")

summary(pix_df$change)

ci <- cwm_res |>
  dplyr::group_by(scenario) |>
  dplyr::summarise(med = median(delta_unw),
                   lo  = quantile(delta_unw, 0.025),
                   hi  = quantile(delta_unw, 0.975))

ggplot(pix_df, aes(x = change, fill = scenario)) +
  geom_histogram(bins = 60, position = "identity", alpha = 0.4) +
  geom_vline(xintercept = 0, colour = "grey30") +
  geom_vline(data = ci, aes(xintercept = med, colour = scenario), linewidth = 0.8) +
  geom_vline(data = ci, aes(xintercept = lo, colour = scenario), linetype = "dashed") +
  geom_vline(data = ci, aes(xintercept = hi, colour = scenario), linetype = "dashed") +
  labs(x = "Change in community borealness (posterior median)",
       y = "Number of pixels",
       fill = NULL, colour = NULL) +
  theme_minimal()
