# 1.- RF modelling ----------
# Modelling GEDI anmd full models with hex_bin data
##------------------------------------------------------------------------#
## @date  2026-05-20
## @update 2026-07-20
## @project c:\\github\\chp_2_structuralForest_resilience
## @R version [Not] Part in a Rumble 4.5.2
## @OS system Windows 10 x64
## @author Jose A. Lastra
## @email jose.lastramunoz@wur.nl | jose.lastra@pucv.cl
##------------------------------------------------------------------------#
# 2.- Libraries --------
pacman::p_load(
  tidyverse,
  sf,
  terra,
  arrow,
  geoarrow,
  spatialRF,
  ranger
)
## Clean environment ----
rm(list = ls(all = T))
gc()
##------------------------------------------------------------------------#
# 3.- Setup  --------
## initial seed -----
SEED <- 42

## outpath -----
outpath <- 'output_models/'
if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

## cluster setup -----
perc <- 0.6
cluster <- floor(parallel::detectCores() * perc) # defining cluster number

## distance threshold ----
distance_thresholds <- c(1000, 5000, 7500, 10000, 12500, 15000, 20000, 50000)

##------------------------------------------------------------------------#
# 4.- Load and data wrangling  --------
## input model data ------
inpath <- 'RF_data/RF_input_data_consolidate.parquet'
data_raw <- inpath %>% open_dataset() %>% st_as_sf()

## data wrangling for modeling -----
### factors to numbers and add FT ------
model_data <- data_raw %>%
  mutate(
    season_num = case_when(
      GS == "2019-20" ~ 1L,
      GS == "2020-21" ~ 2L,
      GS == "2021-22" ~ 3L,
      GS == "2022-23" ~ 4L,
      TRUE ~ NA_real_
    ),
    FT = as.numeric(majority), # identify hex according to mayority
    FT_f = factor(majority, levels = c(1, 2),
                  labels = c('Deciduous', 'Sclerophyllous'))
  ) %>%
  filter(!is.na(GS))

### XY data extraction -----
hex_xy <- model_data %>% st_centroid() %>% st_coordinates()
model_data <- model_data %>% bind_cols(hex_xy)

cat("Available shots:", nrow(model_data), "\n")
cat("Season distribution:\n")
print(table(model_data$GS, model_data$FT_f))

##------------------------------------------------------------------------#
# 5.- Predictor sets  --------
## climate -----
climate_features <- c('pr_zacum', # GS cumulative z-score anomaly precipitation
                      'sm_zacum', # GS cumulative z-score anomaly soil moisture
                      'tmmn_zacum', # GS cumulative z-score anomaly min temp
                      'tmmx_zacum', # GS cumulative z-score anomaly max temp
                      'pdsi_mean' # GS mean pdsi
)

## Stationary (soil - terrain) ------
terrain_soil_features <- c("elev_mean", #  (m)
                           "slope_deg_mean", # slope (°),
                           'sin_aspect_mean', # East-West componente (+ to -)
                           'cos_aspect_mean', # North-South componente (+ to -)
                           'tpi_mean', # Topographic position index
                           'tri_mean', # terrain ruggedness index
                           # soil
                           'bulkd_0_5cm', # bulk density
                           'clay_0_5cm', # clay content
                           'total_awc' # Total Available Water Capacity
)

## GEDI strructure ------
gedi_str <- c("pai_mean", # Plant Area Index total
              "rh_98_mean", # RH 95%
              "agbd_mean", # Biomass (Mg/ha) —
              "cover_mean", # Cover
              "fhd_normal_mean", # Foliage Height Diversity
              paste0('pavd_z_', 1:8, '_mean') # pavd  (mˆ2/mˆ3) between 0 40 (5 m bins)
)

## Global ----
general <- general <- c(
  "FT", # 1 = Deciduous, 2 = Sclerophyllous
  "season_num", # GS numeric
  'evi_cusm' # stress
)


## setup predictors -----
predictor_sets <- list(
  gedi = c(gedi_str, 'FT', 'season_num'),
  full = c(climate_features, terrain_soil_features, gedi_str, general)
)

model_specs <- list(
  gedi              = list(predictors = predictor_sets$gedi, add_interactions = FALSE),
  full              = list(predictors = predictor_sets$full, add_interactions = FALSE)
)

##------------------------------------------------------------------------#
# 6.- Common dataset for modelling  --------
# Union of all predictors, excluding NA rows.

all_predictors <- union(predictor_sets$gedi, predictor_sets$full)

df_model <- model_data %>%
  select(all_of(all_predictors), evi_anom, X, Y) %>%
  drop_na()

cat("Total shots after NA drop (common across all models):", nrow(df_model), "\n")

coords <- df_model %>% select(x = X, y = Y) %>% st_drop_geometry()
distance_matrix <- as.matrix(dist(coords))

##------------------------------------------------------------------------#
# 7.- RF modelling  --------
## running sequentially ----
fitted_models <- list() # empty to save
importance_list <- list() # empty to save

for (model in names(model_specs)) {
  
  cat("\n==== Fitting model:", model, "====\n")
  
  predictors <- model_specs[[model]]$predictors

  ### Data selection ----
  X <- df_model %>% select(evi_anom, all_of(predictors)) %>% st_drop_geometry()
  Y <- df_model$evi_anom
  
  ### model fitting ------
  m_rf <- spatialRF::rf(
    data = X,
    dependent.variable.name = "evi_anom",
    predictor.variable.names = predictors,
    n.cores = cluster,
    distance.matrix = distance_matrix,
    distance.thresholds = distance_thresholds,
    seed = SEED, # change at the beginning of the script
    verbose = TRUE
  )

  ## contribution for model transferibility --
  m_rf <- spatialRF::rf_importance(
    model = m_rf,
    xy = coords, # needs coordinates for cross-validation
    n.cores = cluster, 
    seed = SEED # change at the beginning of the script
  )
  
  ## Export results -----
  ## this generates large size objects
  ## RUN only if needed
  outpath <- 'output_models/'

  outname <- paste0(outpath, 'RF_outputs_', model, '_hex_grid.RData')
  
  save(m_rf,
       file = outname)
  cat("Saved:", outname, "\n")
  
  fitted_models[[model]]  <- m_rf

}

## Cross-model importance table -----
importance_all <- bind_rows(importance_list)
print(importance_all)

## Cross model performance comparison ----
model_comparison <- spatialRF::rf_compare(
  models = fitted_models,
  xy = coords,
  repetitions = 30,
  metrics = c("r.squared", "rmse"),
  seed = SEED,
  n.cores = cluster,
  verbose = TRUE
)

## this generates large size objects
# ## RUN only if needed

outpath <- 'output_models/'
save(fitted_models, model_comparison,
     file = paste0(outpath, 'RF_model_comparison.RData'))










