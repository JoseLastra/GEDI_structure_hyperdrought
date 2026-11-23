# 1.- Cliffs delta analysis ----------
## Cliff delta for the full distribution and by latitude
##------------------------------------------------------------------------#
## @date  2026-03-12
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
  arrow,
  patchwork,
  terra,
  geoarrow,
  effsize,
  broom,
  scales,
  showtext,
  legendry
)
## Clean environment ----
rm(list = ls(all = T))
gc()

##------------------------------------------------------------------------#
# 3.- Input GEDI forest shots -------
## Reading data ------
inpath <- 'GEDI_shots/' 

gedi <- open_dataset(inpath) %>% collect()

## filtering and adding data -----
### Variable selection -----
sel_metrics <- c(
  "shot_number",
  "beam_type",
  "sensitivity",
  "degrade_flag",
  "digital_elevation_model",
  "digital_elevation_model_srtm",
  "elev_lowestmode",
  "cover",
  "fhd_normal",
  "pai",
  "agbd",
  "time",
  'rh_98',
  'rh_50',
  "year",
  "month",
  "ym",
  "DN",
  "EVI_anom",
  "RFD",
  "elev",
  "slope_deg",
  'GS',
  'FT',
  "EVI_state",
  "x",
  "y",
  paste0('pavd_z_', 1:8)
)

### Filtering shots -------
gedi_data <- gedi %>%
  filter(
    !is.na(EVI_anom) &
      rh_98 <= 40 & # excluding potential outliers
      rh_98 >= 5 &
      agbd <= 400 & # excluding potential outliers
      !if_any(
        starts_with(c("pai", "pavd")),
        ~ . < 0
      ) &
      degrade_flag == 0 & # Only good orbits
      sensitivity >= 0.95 &
      month %in% c(10:12, 1:4)
  ) %>%
  mutate(
    FT = factor(
      # Defining forest types
      DN,
      levels = c(1, 2),
      labels = c('Deciduous', 'Sclerophyllous')
    ),
    EVI_state = case_when(
       # 1 affected, 0 unaffected, 3 no data
      EVI_anom < 0 & RFD >= 75 ~ 1L, # Affected forests
      is.na(EVI_anom) & is.na(RFD) ~ 3L, # treating potential NAs
      .default = 0L # Non-affected (default)
    ),
    # defining GS
    GS_start = if_else(month >= 7, year, year - 1),
    GS = paste0(GS_start, "-", substr(GS_start + 1, 3, 4))
  ) %>%
  select(all_of(sel_metrics))

### rename pavd -----
new_names <- paste0("pavd_", seq(0, 35, by = 5), "_", seq(5, 40, by = 5), "_m")
old_names <- paste0('pavd_z_', 1:8)

gedi_data <- gedi_data %>%
  rename_with(~new_names, all_of(old_names))

## Data normalization -----
### Quantile normalization function -----
normalize_p <- function(x) {
  p1 <- quantile(x, 0.01, na.rm = TRUE)
  p99 <- quantile(x, 0.99, na.rm = TRUE)
  if ((p99 - p1) == 0) {
    return(rep(0, length(x)))
  }
  (x - p1) / (p99 - p1)
}

### Normalization -----
gedi_data_norm <- gedi_data %>%
  group_by(FT) %>%
  mutate(across(
    starts_with('rh') |
      starts_with('pai') |
      starts_with('pavd') |
      starts_with('cover') |
      starts_with('fhd') |
      starts_with('agbd'),
    normalize_p,
    .names = '{.col}_norm'
  )) %>%
  ungroup()

## Classify states & filtering spring-summer -----
gedi_filtered <- gedi_data_norm %>%
  filter(
    month %in%
      c(10, 11, 12, 1, 2, 3, 4) &
      time >= '2019-07-01' &
      !is.na(EVI_anom)
  ) %>%
  mutate(
    EVI_state = as.integer(EVI_state),
    EVI_state_f = factor(
      EVI_state,
      levels = c(0, 1),
      labels = c("Non-affected", "Affected")
    ),
    GS = factor(GS)
  )


### Adding GS labels ----
gedi_filtered <- gedi_filtered %>%
  mutate(
    GS_type = case_when(
      GS %in% c("2019-20", "2021-22") ~ "High_stress",
      GS == "2020-21" ~ "Relaxation",
      GS == "2022-23" ~ "Partial",
      TRUE ~ NA_character_
    ),
    GS_type = factor(
      GS_type,
      levels = c("High_stress", "Relaxation", "Partial")
    )
  )

##------------------------------------------------------------------------#
# 4.- General structure analysis -------
## Filtering just relevant columns ------
gedi_structure <- gedi_filtered %>%
  filter(
    !is.na(EVI_state_f)
  )

## Effect size  -----
metrics <- c(
  "rh_98",
  "rh_50",
  "pai",
  "cover",
  "fhd_normal",
  'agbd',
  paste0("pavd_", seq(0, 35, by = 5), "_", seq(5, 40, by = 5), "_m")
)

### Cliff's delta Full distribution -------
cliff_results <- gedi_structure %>%
  group_by(GS, FT) %>%
  group_modify(
    ~ {
      df <- .x

      map_dfr(metrics, function(v) {
        test <- cliff.delta(
          df[[v]] ~ df$EVI_state_f
        )

        tibble(
          metric = v,
          delta = test$estimate,
          conf_l = test$conf.int[1],
          conf_u = test$conf.int[2]
        )
      })
    }
  ) %>%
  ungroup()

# Adding factors
cliff_results <- cliff_results %>%
  mutate(metric = factor(metric, levels = rev(metrics)))

### Cliff's delta latitudinal binning -----
gedi_structure_sf <- gedi_structure %>%
  st_as_sf(coords = c('x', 'y'), crs = 32719) %>%
  st_transform(4326)

# bins
binning <- 1
lat_1 <- c("32–33°S", "33–34°S", "34–35°S", "35–36°S", "36–37°S", "37–38°S")

labels_cut <- lat_1

#### Binning data ---- 
gedi_structure_sf <- gedi_structure_sf %>%
  mutate(
    lat_band = cut(
      st_coordinates(.)[, 2], # lat extraction
      breaks = seq(-38, -32, by = binning), 
      labels = labels_cut %>% rev(),
      include.lowest = TRUE
    )
  )

#### Cliff's delta by bin ------
cliff_results_lat <- gedi_structure_sf %>%
  st_drop_geometry() %>% # no necesitas geometría para el cálculo
  group_by(GS, FT, lat_band) %>%
  group_modify(
    ~ {
      df <- .x

      map_dfr(metrics, function(v) {
        test <- cliff.delta(
          df[[v]] ~ df$EVI_state_f
        )

        tibble(
          metric = v,
          delta = test$estimate,
          conf_l = test$conf.int[1],
          conf_u = test$conf.int[2]
        )
      })
    }
  ) %>%
  ungroup()
