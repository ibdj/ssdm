#### packages ##################################################################

library(tidyverse)
library(myClim)

#### importing tms index (meta)data ##########################################################

# it is only summer temp! 
path <- "~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/tms_data"

path2 <- "~/Library/CloudStorage/OneDrive-GrønlandsNaturinstitut/General - BioBasis/04_GEM_Rådata/TMS4_Logger_Data"

tms_index_mp <- read_delim("~/Library/CloudStorage/OneDrive-Aarhusuniversitet/MappingPlants/02 Modelling future changes/data/tms_index.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE) |> 
  rename(serial = tms_serial)

tms_index_biobasis <- read_csv("~/Library/CloudStorage/OneDrive-GrønlandsNaturinstitut/General - BioBasis/04_GEM_Rådata/TMS4_Logger_Data/tms_metadata2025.csv") |> 
  rename(plot = plot_name, serial = serialnumber, tms_placed = date_installed) |> 
  dplyr::select(-masl)

tms_index_all <- bind_rows(tms_index_mp, tms_index_biobasis) |> 
  dplyr::select(plot,serial, tms_placed, lon, lat, wkt_geom)

#### importing actual tms data with mclim#######################################

files <- c(
  list.files(path,  pattern = "^data_.*\\.csv$", recursive = TRUE, full.names = TRUE),
  list.files(path2, pattern = "^data_.*\\.csv$", recursive = TRUE, full.names = TRUE)
)

length(files)   # sanity check on total count

clean_dir <- file.path(tempdir(), "tms_clean")

dir.create(clean_dir, showWarnings = FALSE)

files_clean <- imap_chr(files, \(f, i) {
  serial <- str_extract(basename(f), "\\d{8}")
  out <- file.path(clean_dir, paste0("data_", serial, "_", i, ".csv"))
  lines <- readLines(f)
  lines <- lines[str_detect(lines, "^\\d+;\\d{4}[.-]\\d{2}[.-]\\d{2} \\d{2}:00")]
  lines |>
    str_replace("^(\\d+;[^;]+);-?\\d+;", "\\1;4;") |>
    str_replace("([^;])$", "\\1;") |>
    writeLines(out)
  out
})

saveRDS(files_clean, "data/files_clean.rds")
readRDS("data/files_clean.rds")

tms <- mc_read_files(files_clean, dataformat_name = "TOMST")
saveRDS(tms, "data/tms.rds")
readRDS("data/tms.rds")

#unique(tms_index_all$serial)

tms <- mc_filter(tms, localities = tms_index_all$serial)
setdiff(locs_before, unique(mc_info(tms)$locality_id))   # who got dropped
"95327113" %in% locs_before                               # was 113 ever read in?

mc_info(tms) |> 
  dplyr::summarise(latest = max(end_date), .by = serial_number)

tms <- mc_join(tms)                    # merges the _1/_2 series per serial
tms <- mc_prep_clean(tms)              # checks step regularity, flags gaps

deploy_crop <- tms_index_all |> 
  dplyr::select(serial, tms_placed) |> 
  dplyr::distinct() |> 
  dplyr::mutate(serial = as.character(serial))

for (i in seq_len(nrow(deploy_crop))) {
  tms <- mc_prep_crop(
    tms,
    start      = as.POSIXct(deploy_crop$tms_placed[i], tz = "UTC"),
    localities = deploy_crop$serial[i]
  )
}

tms <- mc_calc_vwc(tms, soiltype = "universal", output_sensor = "VWC_universal")
tms <- mc_calc_snow(tms, sensor = "TMS_T2")

mc_info(tms)

#### changing in tp a df (from the myclimlist) #################################
# per-logger deployment lookup, keyed to match myclim's serial_number
deploy <- tms_index_all |>
  mutate(serial = as.character(serial)) |>
  dplyr::select(serial, tms_placed, plot, wkt_geom)

deploy |> dplyr::count(serial) |> dplyr::filter(n > 1)

tms_long <- mc_reshape_long(tms)

tms_long_filt <- tms_long |>
  mutate(serial = as.character(serial_number)) |>
  filter(lubridate::minute(datetime) == 0) |>   # hourly first — cuts ~75% immediately
  left_join(deploy, by = "serial") |>
  filter(datetime >= tms_placed)  

str(tms_long_filt)
summary(tms_long_filt)

#### metrics caluclation #######################################################

doy_clim_tms_t1 <- tms_long_filt |>
  filter(sensor_name == "TMS_T1") |>
  mutate(doy = lubridate::yday(datetime)) |>
  group_by(locality_id, doy, plot) |>
  summarise(doy_mean = mean(value, na.rm = TRUE), .groups = "drop")

n_distinct(doy_clim_tms_t1$locality_id)   # expect 69
range(doy_clim_tms_t1$doy)                # expect 1, 366
summary(doy_clim_tms_t1$doy_mean)         # plausible soil temps, no wild outliers

tms_metrics <- doy_clim_tms_t1 |>
  group_by(locality_id, plot) |>
  summarise(
    mat        = mean(doy_mean),
    summer_t   = mean(doy_mean[doy >= 152 & doy <= 243]),  # Jun–Aug
    winter_t   = mean(doy_mean[doy <= 59 | doy >= 335]),   # Dec–Feb
    days_above0 = sum(doy_mean > 0),
    gdd        = sum(pmax(doy_mean, 0)),
    fdd        = sum(pmin(doy_mean, 0)),
    .groups = "drop"
  )

# snow cover pr logger #
record_len <- mc_info(tms) |>
  dplyr::filter(sensor_name == "TMS_T2") |>
  dplyr::summarise(
    years = as.numeric(difftime(max(end_date), min(start_date), units = "days")) / 365.25,
    .by = locality_id, plot
  )

snow_summary <- mc_calc_snow_agg(tms, snow_sensor = "snow")

snow_per_year <- snow_summary |>
  dplyr::left_join(record_len, by = "locality_id") |>
  dplyr::mutate(snow_days_yr = snow_days / years) |>
  dplyr::left_join(
    tms_index_all |> dplyr::mutate(serial = as.character(serial)) |> dplyr::select(serial, plot),
    by = c("locality_id" = "serial")
  ) 

snow_summary

tms_metrics <- tms_metrics |>
  dplyr::left_join(
    tms_index_all |> dplyr::mutate(serial = as.character(serial)) |> dplyr::select(locality_id = serial, plot),
    by = "locality_id"
  )

common1 <- intersect(names(tms_metrics), names(snow_per_year))
common1
tms_metrics <- dplyr::left_join(tms_metrics, snow_per_year, by = common1)
##### fixing geometry ##########################################################

tms_all_sf <- tms_index_all |>
  st_as_sf(wkt = "wkt_geom", crs = 32622)

# check of the coordinates

st_bbox(tms_all_sf)
sum(!st_is_valid(tms_all_sf))
nrow(tms_all_sf)

##### saving outputs ###########################################################

saveRDS(tms_all_sf, "data/tms_all_sf.rds")      # 69 loggers, spatial index
saveRDS(tms_long_filt, "data/tms_long.rds")  
