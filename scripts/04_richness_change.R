library(terra)
library(tidyverse)

# ---- load inputs ----
modelable_species <- readRDS("data/modelable_species.rds")

scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")

borealness <- readRDS("data/idx_modelable.rds") |>
  dplyr::select(species, B_std)

# ---- mean change per species per scenario ----
species_change <- species_change |>
  left_join(name_lookup, by = "taxon") |>
  mutate(species_match = coalesce(species, taxon)) |>
  left_join(idx_modelable |> dplyr::select(species, B_std_new = B_std),
            by = c("species_match" = "species")) |>
  mutate(B_std = coalesce(B_std, B_std_new)) |>
  dplyr::select(-species, -species_match, -B_std_new)

# Add label to outliers
ggplot(species_change, aes(x = B_std, y = mean_change, colour = scenario)) +
  geom_point() +
  ggrepel::geom_text_repel(
    data = species_change |> filter(abs(mean_change) > 0.05),
    aes(label = taxon), size = 2.5, colour = "black"
  ) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  facet_wrap(~scenario) +
  labs(x = "Borealness (B_std)", 
       y = "Mean change in occurrence probability",
       title = "Borealness vs response to warming") +
  theme_minimal()

# Correlation per scenario
species_change |>
  group_by(scenario) |>
  summarise(r = cor(B_std, mean_change, use = "complete.obs"),
            p = cor.test(B_std, mean_change)$p.value)

species_change |>
  group_by(taxon, B_std) |>
  summarise(mean_change_ssp585 = mean_change[scenario == "ssp585"]) |>
  arrange(desc(mean_change_ssp585))

borealness |> pull(species) |> sort()
modelable_species |> sort()

# See what's in borealness that's close to your missing species
borealness |> 
  filter(str_detect(species, "Juncus|Loiseleuria|Ledum|Lycopodium|Deschampsia|Calamagrostis"))

name_lookup <- tibble(
  taxon = c("Juncus trifidus", "Lycopodium annotinum", 
            "Deschampsia flexuosa", "Ledum groenlandicum",
            "Loiseleuria procumbens"),
  species = c("Oreojuncus trifidus", "Spinulum annotinum",
              "Avenella flexuosa", "Rhododendron groenlandicum",
              "Kalmia procumbens")
)
