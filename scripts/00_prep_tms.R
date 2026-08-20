#### packages ##################################################################

library(tidyverse)

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

#### importing actual data #####################################################

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

# BioBasis loggers: lon/lat in WGS84 -> transform to UTM
have_lonlat <- tms_all |>
  dplyr::filter(!is.na(lon)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(32622)

# your loggers: WKT already in UTM 32622
have_wkt <- have_wkt |> dplyr::rename(geometry = wkt_geom)
st_geometry(have_wkt) <- "geometry"