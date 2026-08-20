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
  out <- file.path(clean_dir, paste0("data_", serial, "_", i, ".csv"))  # i keeps it unique
  readLines(f) |>
    str_replace("^(\\d+;[^;]+);-?\\d+;", "\\1;4;") |>
    str_replace("([^;])$", "\\1;") |>
    writeLines(out)
  out
})

tms <- mc_read_files(files_clean, dataformat_name = "TOMST")

mc_info(tms) |> 
  dplyr::summarise(latest = max(end_date), .by = serial_number)

tms <- mc_join(tms)                    # merges the _1/_2 series per serial
tms <- mc_prep_clean(tms)              # checks step regularity, flags gaps
tms <- mc_calc_vwc(tms, soiltype = "universal", output_sensor = "VWC_universal")

mc_info(tms)

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

##### fixing geometry ##########################################################

tms_all_sf <- tms_index_all |>
  st_as_sf(wkt = "wkt_geom", crs = 32622)

# check of the coordinates

st_bbox(tms_all_sf)
sum(!st_is_valid(tms_all_sf))
nrow(tms_all_sf)
