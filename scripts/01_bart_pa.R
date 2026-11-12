#### packages ##################################################################
# install once only - comment out after first run
# devtools::install_github('cjcarlson/embarcadero')

library(embarcadero)
library(pROC)
library(dbarts)
library(raster)
library(parallel)
library(doParallel)
library(foreach)
library(tidyverse)
library(terra)
library(sf)

#### TO DO #####################################################################


#### setup #####################################################################

set.seed(42)

pred_stack <- rast("data/pred_stack.tif")
names(pred_stack)

pred_names <- names(pred_stack)

#### prepare plot data #########################################################

species_frequency <- read_rds("data/species_frequency.rds")

# Get modelable species
modelable_species <- read_rds("data/modelable_species.rds")

species_matrix_pa <- read_rds("data/species_matrix_pa.rds")
names(species_matrix_pa)
rs <- readRDS("data/rs.rds")
names(rs)

# Build predictor + PA matrix
pa_matrix <- species_matrix_pa |>
  as.data.frame() |>
  tibble::rownames_to_column("plot") |>
  dplyr::select(plot, dplyr::all_of(modelable_species)) |>
  dplyr::left_join(sf::st_drop_geometry(rs), by = "plot")

# Verify no NAs in predictors
pa_matrix |>
  dplyr::select(all_of(pred_names)) |>
  summarise(across(everything(), ~sum(is.na(.x)))) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_na") |>
  filter(n_na > 0)

# Predictor matrix for BART
x_train <- pa_matrix |>
  dplyr::select(all_of(pred_names)) |>
  as.data.frame()

# Convert to old raster format for embarcadero/dbarts
pred_rast_stack_r <- raster::stack(pred_stack)

# Dataframe for prediction — keep NAs, track complete cases
pred_df_r    <- as.data.frame(pred_rast_stack_r, na.rm = FALSE)
complete_idx <- complete.cases(pred_df_r)

identical(pred_names, names(pred_stack))

#### BART fitting and spatial projection - parallel ############################

cl <- makeCluster(4)
registerDoParallel(cl)
clusterSetRNGStream(cl, 42)
clusterExport(cl, c("pa_matrix", "x_train", "pred_df_r",
                    "complete_idx", "pred_rast_stack_r", "modelable_species"))

foreach(sp = modelable_species,
        .packages = c("dbarts", "raster"),
        .errorhandling = "pass") %dopar% {
          
          train_data <- cbind(x_train, y = pa_matrix[[sp]]) |> as.data.frame()
          
          model <- dbarts::bart2(
            y ~ .,
            data = train_data,
            keepTrees = TRUE,
            seed = which(modelable_species == sp)
          )
          
          invisible(model$fit$state)   # forces C++ tree state into the R object
          saveRDS(model, paste0("data/bart_pa_", gsub(" ", "_", sp), ".rds"))
          
          # Predict only on complete cases
          pred_vals <- dbarts:::predict.bart(model, 
                                             newdata = pred_df_r[complete_idx, ])
          prob_vals <- colMeans(pred_vals)
          
          # Insert predictions into full vector with NAs for masked pixels
          prob_full <- rep(NA_real_, nrow(pred_df_r))
          prob_full[complete_idx] <- prob_vals
          
          sp_rast_r <- raster::raster(pred_rast_stack_r[[1]])
          raster::values(sp_rast_r) <- prob_full
          
          filename <- paste0("data/sdm_", gsub(" ", "_", sp), ".tif")
          raster::writeRaster(sp_rast_r, filename, overwrite = TRUE)
          
          filename
        }

stopCluster(cl)

# Load saved rasters back into list
species_rasts <- lapply(modelable_species, function(sp) {
  rast(paste0("data/sdm_", gsub(" ", "_", sp), ".tif"))
})
names(species_rasts) <- modelable_species

# Verify they loaded correctly
plot(trim(species_rasts[["Betula nana"]]))

##### 5 fold validation ########################################################

cl <- makeCluster(4)
registerDoParallel(cl)
clusterSetRNGStream(cl, 42)
clusterExport(cl, c("pa_matrix", "x_train", "modelable_species"))

cv_results <- foreach(sp = modelable_species,
                      .packages = c("dbarts", "pROC"),
                      .combine = rbind,
                      .errorhandling = "pass") %dopar% {
                        
                        y <- pa_matrix[[sp]]
                        
                        # Stratified folds
                        pres_idx <- which(y == 1)
                        abs_idx  <- which(y == 0)
                        folds <- numeric(100)
                        folds[pres_idx] <- sample(rep(1:5, length.out = length(pres_idx)))
                        folds[abs_idx]  <- sample(rep(1:5, length.out = length(abs_idx)))
                        
                        fold_auc <- sapply(1:5, function(k) {
                          train_idx <- which(folds != k)
                          test_idx  <- which(folds == k)
                          
                          train_data <- cbind(x_train[train_idx, ], y = y[train_idx]) |> 
                            as.data.frame()
                          
                          cv_model <- dbarts::bart2(
                            y ~ .,
                            data = train_data,
                            keepTrees = TRUE,
                            seed = which(modelable_species == sp) * k
                          )
                          
                          pred <- colMeans(
                            dbarts:::predict.bart(cv_model, 
                                                  newdata = as.data.frame(x_train[test_idx, ]))
                          )
                          
                          as.numeric(auc(roc(y[test_idx], pred, quiet = TRUE)))
                        })
                        
                        data.frame(species = sp, cv_auc = mean(fold_auc), sd_auc = sd(fold_auc))
                      }

stopCluster(cl)

cv_results <- cv_results |>
  arrange(desc(cv_auc)) |>
  mutate(assessment = case_when(
    cv_auc >= 0.90 ~ "Excellent",
    cv_auc >= 0.80 ~ "Good",
    cv_auc >= 0.70 ~ "Acceptable",
    cv_auc >= 0.60 ~ "Poor",
    TRUE           ~ "Fail"
  ))

cv_results |> 
  tibble::as_tibble() |> 
  print(n = Inf)


#### learning curves ###########################################################
sizes <- c(30, 50, 70, 90, 110, 130)

cl <- makeCluster(4)
registerDoParallel(cl)
clusterSetRNGStream(cl, 42)
clusterExport(cl, c("pa_matrix", "x_train", "modelable_species"))

learn_curves <- foreach(sp = modelable_species,
                        .packages = c("dbarts", "pROC"),
                        .combine = rbind,
                        .errorhandling = "remove") %dopar% {
                          
                          y <- pa_matrix[[sp]]
                          
                          do.call(rbind, lapply(sizes, function(n) {
                            aucs <- sapply(1:5, function(rep) {
                              idx_train <- sample(seq_along(y), n)
                              idx_test  <- setdiff(seq_along(y), idx_train)
                              if (length(unique(y[idx_train])) < 2 ||
                                  length(unique(y[idx_test]))  < 2) return(NA)
                              
                              m <- dbarts::bart2(y ~ .,
                                                 data = cbind(x_train[idx_train, ], y = y[idx_train]) |> as.data.frame(),
                                                 keepTrees = TRUE, seed = rep)
                              p <- colMeans(dbarts:::predict.bart(m,
                                                                  newdata = as.data.frame(x_train[idx_test, ])))
                              as.numeric(auc(roc(y[idx_test], p, quiet = TRUE)))
                            })
                            data.frame(species = sp, n_train = n, mean_auc = mean(aucs, na.rm = TRUE))
                          }))
                        }
stopCluster(cl)

ggplot(learn_curves, aes(n_train, mean_auc, colour = species)) +
  geom_line() +
  labs(x = "Training sample size", y = "Mean AUC", title = "Learning curves") +
  theme_minimal()

#### diagnostics ###############################################################
# Load the fitted models from the parallel loop (same models as the rasters)
bart_models <- setNames(lapply(modelable_species, function(sp) {
  train_data <- cbind(x_train, y = pa_matrix[[sp]]) |> as.data.frame()
  dbarts::bart2(y ~ ., data = train_data, keepTrees = TRUE,
                seed = which(modelable_species == sp))
}), modelable_species)

saveRDS(bart_models, "data/bart_models.rds")
length(bart_models[["Juncus trifidus"]]$fit$state)   # should be 4

# Training AUC
auc_results <- map_dfr(modelable_species, function(sp) {
  predicted <- colMeans(
    dbarts:::predict.bart(bart_models[[sp]], 
                          newdata = as.data.frame(x_train))
  )
  observed <- pa_matrix[[sp]]
  
  auc_val <- auc(roc(observed, predicted, quiet = TRUE))
  
  tibble(
    species = sp,
    n_plots = sum(observed),
    auc = round(as.numeric(auc_val), 3)
  )
}) |>
  arrange(desc(auc))


invisible(model$fit$state)
saveRDS(model, paste0("data/bart_pa_", gsub(" ", "_", sp), ".rds"))


## per-draw predictions

sp <- "Juncus trifidus"   # adjust to your naming
model <- bart_models[[sp]]

pred_env <- as.data.frame(na.omit(terra::values(pred_stack)))

draws_cur <- dbarts:::predict.bart(model, newdata = pred_env)
dim(draws_cur)   # should be n_draws x n_cells (e.g. 2000 x ~cells)

# per draw mean one sceneario at a time: 

delta <- deltas[["ssp585"]]   # adjust to how your deltas are stored

pred_ssp <- pred_env
pred_ssp$summer_t_int <- pred_ssp$summer_t_int + delta

draws_ssp <- dbarts:::predict.bart(model, newdata = pred_ssp)

change_draws <- rowMeans(draws_ssp) - rowMeans(draws_cur)   # mean change per draw
mean(change_draws)
quantile(change_draws, c(0.025, 0.975))

# table of all credible intervals for alle species

cred_table <- map_dfr(modelable_species, function(sp) {
  model <- bart_models[[sp]]
  draws_cur <- dbarts:::predict.bart(model, newdata = pred_env)
  cur <- rowMeans(draws_cur)
  
  map_dfr(names(deltas), function(scen) {
    pred_ssp <- pred_env
    pred_ssp$summer_t_int <- pred_ssp$summer_t_int + deltas[[scen]]
    change <- rowMeans(dbarts:::predict.bart(model, newdata = pred_ssp)) - cur
    
    tibble(species = sp, scenario = scen,
           mean_change = mean(change),
           lower95 = quantile(change, 0.025),
           upper95 = quantile(change, 0.975),
           p_increase = mean(change > 0))
  })
})
print(cred_table, n = Inf)
 #### explanatory ###############################################################

#### variable importance #######################################################

varimp_long <- map_dfr(modelable_species, function(sp) {
  vc <- bart_models[[sp]]$varcount
  
  # Average across chains (dim 1) and iterations (dim 2)
  counts <- apply(vc, 3, mean)
  
  tibble(
    variable = pred_names,
    importance = counts / sum(counts),
    species = sp
  )
})

head(varimp_long)

ggplot(varimp_long, aes(x = variable, y = species, fill = importance)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Importance") +
  labs(x = "Variable", y = "Species",
       title = "Variable importance across BART SDMs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"))

#### VERSION TWO ##################################################################

# Calculate full model R2 per species
r2_results <- map_dfr(modelable_species, function(sp) {
  predicted <- colMeans(
    dbarts:::predict.bart(bart_models[[sp]], 
                          newdata = as.data.frame(x_train))
  )
  observed <- pa_matrix[[sp]]
  
  ss_res <- sum((observed - predicted)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  r2 <- 1 - (ss_res / ss_tot)
  
  tibble(species = sp, r2 = round(r2, 3))
}) |>
  arrange(desc(r2))

# R2 drop when each variable is excluded
varimp_r2 <- map_dfr(modelable_species, function(sp) {
  observed <- pa_matrix[[sp]]
  
  # Full model R2
  full_pred <- colMeans(dbarts:::predict.bart(bart_models[[sp]], 
                                              newdata = as.data.frame(x_train)))
  ss_tot <- sum((observed - mean(observed))^2)
  r2_full <- 1 - sum((observed - full_pred)^2) / ss_tot
  
  # R2 drop per variable
  map_dfr(pred_names, function(var) {
    x_reduced <- x_train
    x_reduced[[var]] <- mean(x_train[[var]], na.rm = TRUE)
    
    red_pred <- colMeans(dbarts:::predict.bart(bart_models[[sp]], 
                                               newdata = as.data.frame(x_reduced)))
    r2_reduced <- 1 - sum((observed - red_pred)^2) / ss_tot
    
    tibble(species = sp, variable = var, 
           r2_drop = max(0, r2_full - r2_reduced))
  })
})

# Normalise r2_drop to sum to r2_full per species
varimp_r2_norm <- varimp_r2 |>
  left_join(r2_results, by = "species") |>
  group_by(species) |>
  mutate(r2_drop_norm = (r2_drop / sum(r2_drop)) * r2) |>
  ungroup()

# Recalculate unexplained
unexplained <- r2_results |>
  mutate(variable = "Unexplained",
         r2_drop_norm = 1 - r2)

# Combine
plot_data <- varimp_r2_norm |>
  dplyr::select(species, variable, r2_drop_norm) |>
  bind_rows(dplyr::select(unexplained, species, variable, r2_drop_norm)) |>
  mutate(
    species = factor(species, levels = r2_results$species),
    variable = factor(variable, levels = c("Unexplained", rev(pred_names)))
  )

# Plot
ggplot(plot_data, aes(x = r2_drop_norm, y = species, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(
    values = c(
      "elevation"      = "#2ca02c",
      "slope"          = "#ff7f0e",
      "ndwi"           = "#1f77b4",
      "hli"            = "#9467bd",
      "snowfree"       = "#e377c2",
      "summer_t_int"   = "#d62728",
      "Unexplained"    = "#D9D9D9"
    )
  ) +
  labs(x = "Proportion of variance", y = NULL,
       fill = "Variable",
       title = "Explained variance partitioned by predictor") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))

##### variables ################################################################

# Pivot x_train to long format for faceted plotting
x_train_long <- x_train |>
  pivot_longer(everything(), 
               names_to = "variable", 
               values_to = "value")

# Histogram of plot distribution across each predictor
ggplot(x_train_long, aes(x = value)) +
  geom_histogram(bins = 20, fill = "#2E75B6", colour = "white") +
  facet_wrap(~ variable, scales = "free") +
  labs(x = "Value", y = "Number of plots",
       title = "Distribution of plots across environmental predictors") +
  theme_minimal()


#plotting the table

cred_table <- cred_table |>
  mutate(species = fct_reorder(species, mean_change, .fun = max, .desc = TRUE))

ggplot(cred_table, aes(x = species, y = mean_change,
                       colour = scenario, shape = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = lower95, ymax = upper95),
                  position = position_dodge(width = 0.7), size = 0.3) +
  labs(x = NULL, y = "Change in mean occurrence probability",
       colour = "Scenario", shape = "Scenario") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
