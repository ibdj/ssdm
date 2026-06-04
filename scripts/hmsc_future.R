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
library(ggplot2)
library(corrplot)
library(viridis)

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

# X matrix - environmental predictors (without ndvi)
X <- abiotic_plot |>
  dplyr::select(elevation, slope, aspect_sin, aspect_cos,
                twi, temp_predicted) |>
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
  twi + temp_predicted

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
                nParallel = nChains,   # <-- to run in parallel
                verbose   = 500
)

# Save model immediately after sampling
saveRDS(m, "data/hmsc_model_pa.rds")

summary(m)


# ==============================================================================

# ==============================================================================
# HMSC Model Visualizations
# Requires: Hmsc, ggplot2, corrplot, tidyr, dplyr, viridis
# ============================================================

# ---- Load your saved model ----
m <- readRDS("data/hmsc_model_pa.rds")

# ==============================================================
# 1. MCMC CONVERGENCE DIAGNOSTICS
# Check that chains have converged before interpreting results
# ==============================================================

# Compute posterior samples
mpost <- convertToCodaObject(m)

# Potential Scale Reduction Factor (Gelman-Rubin)
# Values < 1.1 indicate good convergence
psrf_beta <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
psrf_gamma <- gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf

cat("=== MCMC Convergence (Gelman-Rubin PSRF) ===\n")
cat("Beta parameters - Mean PSRF:", round(mean(psrf_beta[,1]), 3), "\n")
cat("Beta parameters - Max PSRF:", round(max(psrf_beta[,1]), 3), "\n")
cat("Gamma parameters - Mean PSRF:", round(mean(psrf_gamma[,1]), 3), "\n")

# Effective sample sizes
ess_beta <- effectiveSize(mpost$Beta)
cat("Beta ESS - Mean:", round(mean(ess_beta), 1), "\n\n")

# PSRF Heatmap
psrf_df <- data.frame(
  parameter = names(psrf_beta[,1]),
  psrf = psrf_beta[,1]
)

png("hmsc_1_convergence.png", width = 900, height = 500, res = 120)
ggplot(psrf_df, aes(x = reorder(parameter, psrf), y = psrf, fill = psrf)) +
  geom_col() +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_gradient2(low = "#2196F3", mid = "#FFF176", high = "#F44336",
                       midpoint = 1.05, limits = c(1, max(psrf_df$psrf))) +
  coord_flip() +
  labs(title = "MCMC Convergence — Gelman-Rubin PSRF",
       subtitle = "Dashed line at 1.1 (values below = good convergence)",
       x = NULL, y = "PSRF", fill = "PSRF") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")
dev.off()
cat("Saved: hmsc_1_convergence.png\n")


# ==============================================================
# 2. BETA COEFFICIENTS — Species responses to environment
# Shows which species respond positively/negatively to each covariate
# ==============================================================

# Compute posterior means and support (probability of positive effect)
postBeta <- getPostEstimate(m, parName = "Beta")
beta_mean <- postBeta$mean
beta_support <- postBeta$support  # P(effect > 0)

# Plot: Beta means as heatmap
rownames(beta_mean) <- m$covNames
colnames(beta_mean) <- m$spNames

beta_df <- as.data.frame(beta_mean) %>%
  rownames_to_column("covariate") %>%
  pivot_longer(-covariate, names_to = "species", values_to = "estimate")

png("hmsc_2_beta_heatmap.png", width = 1000, height = 600, res = 120)
ggplot(beta_df, aes(x = species, y = covariate, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#B71C1C",
                       midpoint = 0, name = "β estimate") +
  labs(title = "Species–Environment Relationships (Beta coefficients)",
       subtitle = "Red = positive response, Blue = negative response",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank())
dev.off()
cat("Saved: hmsc_2_beta_heatmap.png\n")


# ==============================================================
# 3. BETA SUPPORT — Statistical support for effects
# Shows which coefficients have strong evidence (>0.9 or <0.1)
# ==============================================================

rownames(beta_support) <- m$covNames
colnames(beta_support) <- m$spNames

support_df <- as.data.frame(beta_support) %>%
  rownames_to_column("covariate") %>%
  pivot_longer(-covariate, names_to = "species", values_to = "support") %>%
  mutate(significant = case_when(
    support > 0.9 ~ "Positive (>0.9)",
    support < 0.1 ~ "Negative (<0.1)",
    TRUE ~ "Uncertain"
  ))

png("hmsc_3_beta_support.png", width = 1000, height = 600, res = 120)
ggplot(support_df, aes(x = species, y = covariate, fill = significant)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Positive (>0.9)" = "#C62828",
                               "Negative (<0.1)" = "#1565C0",
                               "Uncertain" = "grey90"),
                    name = "Support") +
  labs(title = "Statistical Support for Species–Environment Effects",
       subtitle = "Strong positive (P>0.9) or negative (P<0.1) effects highlighted",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid = element_blank())
dev.off()
cat("Saved: hmsc_3_beta_support.png\n")


# ==============================================================
# 4. VARIANCE PARTITIONING
# How much of species occurrence is explained by each predictor?
# ==============================================================

VP <- computeVariancePartitioning(m)

# Plot: stacked bar chart of R2 contributions per species
vp_df <- as.data.frame(VP$vals)
colnames(vp_df) <- m$spNames
vp_df$predictor <- rownames(vp_df)
vp_long <- pivot_longer(vp_df, -predictor, names_to = "species", values_to = "r2")

# Average VP across species for summary
vp_avg <- VP$vals %>%
  rowMeans() %>%
  data.frame(predictor = rownames(VP$vals), avg_r2 = .)

png("hmsc_4_variance_partitioning.png", width = 1100, height = 600, res = 120)
ggplot(vp_long, aes(x = species, y = r2, fill = predictor)) +
  geom_col(position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_viridis_d(option = "turbo", name = "Predictor") +
  labs(title = "Variance Partitioning",
       subtitle = "Proportion of species occurrence explained by each predictor",
       x = NULL, y = "R²") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "right")
dev.off()
cat("Saved: hmsc_4_variance_partitioning.png\n")


# ==============================================================
# 5. SPECIES CO-OCCURRENCE (OMEGA matrix)
# Residual correlations between species after accounting for environment
# Positive = co-occur more than expected; Negative = avoid each other
# ==============================================================

OmegaCor <- computeAssociations(m)
supportLevel <- 0.95

# Extract mean omega and support
toPlot <- ((OmegaCor[[1]]$support > supportLevel) +
             (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) *
  OmegaCor[[1]]$mean

rownames(toPlot) <- m$spNames
colnames(toPlot) <- m$spNames

png("hmsc_5_omega_correlations.png", width = 900, height = 800, res = 120)
corrplot(toPlot,
         method = "color",
         col = colorRampPalette(c("#1565C0", "white", "#B71C1C"))(200),
         tl.cex = 0.85,
         tl.col = "black",
         title = paste0("Residual Species Associations (>", supportLevel*100, "% support)"),
         mar = c(0, 0, 2, 0),
         addCoef.col = NULL,
         diag = FALSE)
dev.off()
cat("Saved: hmsc_5_omega_correlations.png\n")


# ==============================================================
# 6. GAMMA COEFFICIENTS — Trait effects on species responses
# How do species traits modify their responses to environment?
# ==============================================================

postGamma <- getPostEstimate(m, parName = "Gamma")
gamma_mean    <- postGamma$mean
gamma_support <- postGamma$support

# Gamma is structured as [covariates × traits] in Hmsc
rownames(gamma_mean)    <- m$covNames   # was m$trNames — wrong
colnames(gamma_mean)    <- m$trNames    # was m$covNames — wrong
rownames(gamma_support) <- m$covNames
colnames(gamma_support) <- m$trNames

gamma_df <- as.data.frame(gamma_mean) %>%
  rownames_to_column("trait") %>%
  pivot_longer(-trait, names_to = "covariate", values_to = "estimate")

gamma_support_df <- as.data.frame(gamma_support) %>%
  rownames_to_column("trait") %>%
  pivot_longer(-trait, names_to = "covariate", values_to = "support")

gamma_df$support <- gamma_support_df$support
gamma_df$sig <- ifelse(gamma_df$support > 0.9 | gamma_df$support < 0.1, "*", "")

png("hmsc_6_gamma_trait_env.png", width = 900, height = 500, res = 120)
ggplot(gamma_df, aes(x = trait, y = covariate, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 5, color = "black") +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#B71C1C",
                       midpoint = 0, name = "γ estimate") +
  labs(title = "Trait–Environment Relationships (Gamma coefficients)",
       subtitle = "* = >90% posterior support | Red = positive, Blue = negative",
       x = "Covariate", y = "Trait") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
cat("Saved: hmsc_6_gamma_trait_env.png\n")


# ==============================================================
# 7. PREDICTIVE PERFORMANCE
# Cross-validated model fit statistics (AUC for P/A data)
# ==============================================================

# Note: this is computationally intensive — reduce k for large models
# Use partition if already computed, otherwise compute here
colnames(m$studyDesign)
partition <- createPartition(m, nfolds = 2, column = "plot")
preds <- computePredictedValues(m, partition = partition)

MF <- evaluateModelFit(hM = m, predY = preds)

perf_df <- data.frame(
  species  = m$spNames,
  AUC      = MF$AUC,
  TjurR2   = MF$TjurR2,
  RMSE     = MF$RMSE
)

cat("\n=== Predictive Performance ===\n")
print(perf_df)

png("hmsc_7_model_fit_auc.png", width = 900, height = 500, res = 120)
ggplot(perf_df, aes(x = reorder(species, AUC), y = AUC, fill = AUC)) +
  geom_col() +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "#616161") +
  scale_fill_gradient(low = "#FFCCBC", high = "#B71C1C") +
  coord_flip() +
  labs(title = "Cross-validated AUC per Species",
       subtitle = "Dashed line at AUC = 0.7 (acceptable discrimination)",
       x = NULL, y = "AUC (2-fold CV)", fill = "AUC") +
  theme_minimal(base_size = 11)
dev.off()
cat("Saved: hmsc_7_model_fit_auc.png\n")


# ==============================================================
# 8. RESPONSE CURVES — Marginal species responses to covariates
# Predicted occurrence probability across a gradient of one variable
# (holding all others at their mean)
# ==============================================================

# Adjust "focal_covariate" to the name of your variable of interest
focal_covariate <- m$covNames[2]  # Change index as needed

Gradient <- constructGradient(m,
                              focalVariable = focal_covariate,
                              non.focalVariables = list(type = "mean"))

predY <- predict(m,
                 XData = Gradient$XDataNew,
                 studyDesign = Gradient$studyDesignNew,
                 ranLevels = Gradient$rLNew,
                 expected = TRUE)

plotGradient(m, Gradient, predY, measure = "Y",
             showData = TRUE, jigger = 0.1,
             las = 2)

# Save response curve plot
png("hmsc_8_response_curves.png", width = 1100, height = 700, res = 120)
plotGradient(m, Gradient, predY, measure = "Y",
             showData = TRUE, jigger = 0.1, las = 2,
             main = paste("Species response curves along", focal_covariate))
dev.off()
cat("Saved: hmsc_8_response_curves.png\n")


# ==============================================================
# SUMMARY
# ==============================================================
cat("\n=== All plots saved ===\n")
cat("1. hmsc_1_convergence.png         — MCMC diagnostics\n")
cat("2. hmsc_2_beta_heatmap.png        — Beta coefficients\n")
cat("3. hmsc_3_beta_support.png        — Beta statistical support\n")
cat("4. hmsc_4_variance_partitioning.png — Variance explained per predictor\n")
cat("5. hmsc_5_omega_correlations.png  — Residual species co-occurrence\n")
cat("6. hmsc_6_gamma_trait_env.png     — Trait x environment (Gamma)\n")
cat("7. hmsc_7_model_fit_auc.png       — Cross-validated AUC\n")
cat("8. hmsc_8_response_curves.png     — Species response curves\n")

