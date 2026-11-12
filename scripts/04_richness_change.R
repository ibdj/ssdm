library(terra)
library(tidyverse)

# ---- load inputs ----
modelable_species <- readRDS("data/modelable_species.rds")

scenarios <- c("current", "ssp126", "ssp245", "ssp370", "ssp585")

borealness <- readRDS("data/idx_modelable.rds") |>
  dplyr::select(species, B_std)

# ---- mean change per species per scenario ----
species_change <- readRDS("data/boreal_sp_change_future.rds")

names(species_change)
sum(is.na(species_change$B_std))   # 0 expected

# Add label to outliers
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

# richness

richness <- purrr::map(scenarios, \(s) {
  f <- file.path("data", paste0("sdm_", gsub(" ", "_", modelable_species),
                                if (s == "current") "" else paste0("_", s), ".tif"))
  sum(terra::rast(f))
}) |> terra::rast()
names(richness) <- scenarios

rich_change <- richness[[-1]] - richness[["current"]]

terra::global(richness[["current"]], c("min", "max", "mean"), na.rm = TRUE)

tibble::tibble(
  scenario     = scenarios[-1],
  mean_change  = terra::global(rich_change, "mean", na.rm = TRUE)[[1]],
  sd_change    = terra::global(rich_change, "sd",   na.rm = TRUE)[[1]],
  pct_increase = 100 * terra::global(rich_change > 0, "mean", na.rm = TRUE)[[1]],
  pct_decrease = 100 * terra::global(rich_change < 0, "mean", na.rm = TRUE)[[1]]
)

# credibility of the model 

pred_stack        <- terra::rast("data/pred_stack.tif")
modelable_species <- readRDS("data/modelable_species.rds")
deltas <- c(ssp126 = 0.9428791, ssp245 = 1.9663557,
            ssp370 = 2.8462466, ssp585 = 3.8888474)

set.seed(42)
pix <- terra::spatSample(pred_stack, size = 10000, method = "random",
                         na.rm = TRUE, as.df = TRUE)

newdata <- c(list(current = pix),
             purrr::map(deltas, \(d) dplyr::mutate(pix, summer_t_int = summer_t_int + d)))

nrow(pix)                      # 10000 expected
purrr::map_dbl(newdata, \(x) mean(x$summer_t_int))

#step 2 credibility
draws <- purrr::map_dfr(modelable_species, \(sp) {
  m <- readRDS(paste0("data/bart_pa_", gsub(" ", "_", sp), ".rds"))
  purrr::map_dfr(names(newdata), \(s) {
    p <- dbarts:::predict.bart(m, newdata = newdata[[s]])   # draws x pixels
    tibble::tibble(taxon = sp, scenario = s,
                   draw = seq_len(nrow(p)), p_mean = rowMeans(p))
  })
}, .progress = TRUE)

saveRDS(draws, "data/richness_draws_future.rds")

dplyr::count(draws, taxon, scenario) |> dplyr::distinct(n)   # one value expected

# step 3

rich_draws <- draws |>
  dplyr::summarise(richness = sum(p_mean), .by = c(scenario, draw)) |>
  tidyr::pivot_wider(names_from = scenario, values_from = richness) |>
  tidyr::pivot_longer(dplyr::starts_with("ssp"), names_to = "scenario", values_to = "future") |>
  dplyr::mutate(change = future - current)

rich_ci <- rich_draws |>
  dplyr::summarise(mean   = mean(change),
                   lower  = quantile(change, 0.025),
                   upper  = quantile(change, 0.975),
                   p_gain = mean(change > 0),
                   .by = scenario)

saveRDS(rich_ci, "data/richness_ci_future.rds")
rich_ci
