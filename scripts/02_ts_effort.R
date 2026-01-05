################################################################################
# title
################################################################################
#
# Juan Carlos Villaseñor-Derbez
# jc_villasenor@miami.edu
# date
#
# Build a brief figure of before, duirng, after fishing effort using GFW and 
# NOAA hurricane data.
#
################################################################################

## SET UP ######################################################################
################################################################################
# title: Time series of fishing effort relative to hurricane exposure
# author: Juan Carlos Villaseñor-Derbez <jc_villasenor@miami.edu>
# date: (update when running)
#
# Purpose:
# Build a concise figure that compares fishing effort (Global Fishing Watch)
# before, during and after hurricane exposure using NOAA/IBTrACS hurricane tracks
# and a storm-wind exposure model. Output is saved to `results/img/ts_effort.png`.
#
# Notes:
# - Uses `gfwr` to download daily fishing-raster tiles for the Gulf of Mexico.
# - Uses `stormwindmodel` to calculate grid-level wind exposure time series.
# - Effort metrics are normalized by pixel area (km^2) before aggregation.
################################################################################

# Load packages ----------------------------------------------------------------
library(here)
## SET UP ######################################################################

# Load packages ----------------------------------------------------------------
library(gfwr)
library(stormwindmodel)
library(tidyverse)
library(sf)
library(terra)

# Define user functions --------------------------------------------------------
# Function to fix the length of timestamp strings from Hurricane data
# Inputs:
# - char: character vector of numeric components (month, day, hour, minute)
# - desired_length: integer target width (e.g. 2 for '03')
# Returns: character vector padded with leading zeroes to `desired_length`
fix_length <- function(char, desired_length){
  current_length <- str_length(char)
  
  missing <- desired_length - current_length
  
  zeroes <- character(length = length(missing))
  
  zeroes[which(missing >= 1)] <- 0
  
  char <- paste0(zeroes, char)
  
  return(char)
}

# Compute time series of sustained wind speed for each grid cell
# Inputs:
# - track: a single hurricane track (data.frame or sf) for one storm
# - grid: data.frame with grid coordinates and `gridid` used by stormwindmodel
# Returns: tibble with columns `grid_id`, `lon`, `lat`, `date`, `vmax_sust`
get_hurricane_exposure <- function(track, grid) {
  # Calculate exposure using stormwindmodel; this returns a matrix-like
  # object with times as rownames and grid cells as columns.
  exposure <- stormwindmodel::calc_grid_winds(hurr_track = track,
                                              grid_df = grid)
  # Convert the exposure object into a long tibble with one row per
  # (time, grid cell) and keep only positive sustained wind values.
  ts <- exposure[["vmax_sust"]] |> 
    tibble::as_tibble(rownames = "time") |> 
    tidyr::pivot_longer(cols = -time, 
                        names_to = "grid_id",
                        values_to = "vmax_sust") |> 
    dplyr::mutate(time = lubridate::ymd_hms(time),
                  date = lubridate::date(time)) |> 
    dplyr::left_join(grid, by = dplyr::join_by(grid_id == gridid)) |> 
    dplyr::group_by(grid_id, glon, glat, date) |> 
    dplyr::summarize(vmax_sust = max(vmax_sust, na.rm = T),
                     .groups = "drop") |> 
    dplyr::select(grid_id, lon = glon, lat = glat, date, vmax_sust) |> 
    dplyr::filter(vmax_sust > 0)

  # Return the time series tibble
  return(x = ts)
}

# Load spatial mask / study region --------------------------------------------
# `World_Seas_IHO_v3` shapefile is stored locally in Box; we simplify the
# polygon to speed spatial operations. `GoM` is an `sf` polygon for the Gulf
# of Mexico used to filter tracks and region-based downloads.
GoM <- read_sf("/Users/jcvd/Library/CloudStorage/Box-Box/01_project_data/central/World_Seas_IHO_v3") |> 
  filter(NAME == "Gulf of Mexico") |> 
  st_simplify(dTolerance = 10000)

# Hurricanes: read IBTrACS (NOAA) CSV and keep full column names by reading
# a header row separately. This dataset contains best-track positions and wind
# metrics for all storms in the North Atlantic.
col_names <- names(read_csv("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.NA.list.v04r01.csv", n_max = 0))
hur <- read_csv("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.NA.list.v04r01.csv",
                col_names = col_names,
                skip = 2)

## PROCESSING ##################################################################
## Helper: download daily raster tiles for a single year ----------------------
## This wraps `gfwr::get_raster()` to fetch Global Fishing Watch daily rasters
## for the Gulf of Mexico, then performs light renaming and creates a
## `grid_id` identifier used to match pixels to hurricane exposure.
my_get_raster <- function(yr) {
  effort <- get_raster(spatial_resolution = "LOW",
      temporal_resolution = "DAILY",
      region = GoM,
      region_source = "USER_SHAPEFILE",
      start_date = paste0(yr, "-06-01"),
      end_date = paste0(yr, "-12-01"),
      filter_by = "flag = 'USA'") |> 
    rename(lat = Lat,
      lon = Lon,
      date = `Time Range`,
      effort_n_vessels = `Vessel IDs`,
      effort_hours = `Apparent Fishing Hours`) |> 
    # shift pixel centers slightly and build a joinable id
    mutate(lat = lat + 0.05,
      lon = lon + 0.05,
      grid_id = paste(lon, lat))

  return(effort)
}

## Download effort rasters for each year in the analysis window
rast <- map_dfr(2015:2020,
                my_get_raster)

# PROCESSING ###################################################################

## Prepare hurricanes ----------------------------------------------------------
gulf_hur <- hur |> 
  filter(SEASON >= 2015,
         !NAME == "UNNAMED") |> 
  mutate(ID = paste0(SEASON, "_", NAME)) |> 
  select(ID, LON, LAT) |> 
  st_as_sf(coords = c("LON", "LAT"),
           crs = "EPSG:4326") |> 
  group_by(ID) |> 
  summarize(n = n(),
            do_union = F) |> 
  st_cast("LINESTRING") |> 
  st_filter(GoM) |> 
  pull(ID) |> 
  unique()
         
hurricanes <- hur |> 
  select(name = NAME, year = SEASON, date = ISO_TIME,
         latitude = LAT, longitude = LON,
         wind = USA_WIND, sshs = USA_SSHS) |> 
  mutate(ID = paste0(year, "_", name)) |> 
  filter(ID %in% gulf_hur) |> 
  mutate(date = paste0(year(date),
                       fix_length(month(date), desired_length = 2),
                       fix_length(day(date), desired_length = 2),
                       fix_length(hour(date), desired_length = 2),
                       fix_length(minute(date), desired_length = 2)))

## Calculate hurricane exposure for each ocean pixel -------------------------
## Build a unique grid of pixel centers used in the raster download and then
## convert to a simple data.frame matching the `stormwindmodel` grid format.
ocean_grid <- rast |> 
  select(glon = lon, glat = lat) |> 
  distinct() |> 
  mutate(gridid = paste(glon, glat),
         glandsea = F) |> 
  as.data.frame()

## Estimate pixel area (km^2) from a raster of cell sizes. This area is used
## to normalize vessel counts / hours to per-km2 units.
ocean_grid_rast <- rast(ocean_grid |>
                          select(1:2),
                        crs = "EPSG:4326") |> 
  cellSize() |> 
  as.data.frame(xy = T) |> 
  rename(lon = x, lat = y) |> 
  mutate(area = area / 1e6) |> 
  as_tibble()

## Use parallel workers (mirai) to speed up exposure calculations across storms
mirai::daemons(20)
## For each hurricane (split by ID) compute grid-level sustained-wind time
## series using `get_hurricane_exposure()` and combine into one data frame.
exposure <- hurricanes %>%
  split(.$ID) |> 
  map_dfr(in_parallel(\(x) get_hurricane_exposure(x,
                                                  grid = grid),
                      get_hurricane_exposure = get_hurricane_exposure,
                      grid = ocean_grid),
          .id = "hurricane")

sshs_hurricanes <- hurricanes |> 
  group_by(ID) |> 
  summarize(sshs = max(sshs))

panel <- exposure |> 
  filter(vmax_sust > 18) |>
  group_by(hurricane, grid_id, lon, lat) |>
  summarize(first_day = min(date),
            last_day = max(date),
            .groups = "drop") |>
  mutate(start_pre = first_day - 60,
         end_post = last_day + 45) |>
  mutate(date = map2(start_pre, end_post, ~ seq(.x, .y, by = "day"))) |>
  unnest(date) |>
  mutate(period = case_when(between(date, first_day, last_day) ~ "during",
                            date < first_day ~ "pre",
                            date > last_day ~ "post"),
         period = fct_relevel(period, c("pre", "during", "post")),
         rel_time_left = date - first_day,
         rel_time_right = date - last_day,
         rel_time = case_when(period == "pre" ~ as.numeric(rel_time_left),
                              period == "during" ~ 0,
                              period == "post" ~ as.numeric(rel_time_right))) |>
  left_join(exposure, by = join_by(hurricane, grid_id, lon, lat, date)) |>
  left_join(rast |> select(-c(lat, lon)), by = join_by(grid_id, date)) |>
  left_join(ocean_grid_rast, by = join_by(lon, lat)) |> 
  left_join(sshs_hurricanes, by = join_by(hurricane == ID)) |> 
  select(hurricane, grid_id, lon, lat, area, date, rel_time, period, sshs, vmax_sust, effort_n_vessels, effort_hours) |> 
  replace_na(list(v_max_sust = 0,
                  effort_n_vessels = 0,
                  effort_hours = 0)) |> 
  mutate(effort_n_vessels_km2 = effort_n_vessels / area,
         effort_hours_km2 = effort_hours / area)

## Aggregate statistics for plotting -----------------------------------------
## Compute the mean and standard error of vessel counts (per km^2) by
## the three analysis periods ('pre', 'during', 'post'). These values are
## used to draw a reference horizontal line (pre-storm mean) in the plot.
period <- group_by(panel,
         period) |> 
  summarize(mean = mean(effort_n_vessels_km2,
         na.rm = T),
       se = sd(effort_n_vessels_km2, na.rm = T) / sqrt(n()))

## Plotting: visualize mean vessel activity relative to storm exposure
## - x: days relative to first/last storm-force wind day
## - y: vessel count per km^2
## - vertical dashed line marks the storm period boundary (rel_time == 0)
## - horizontal line shows the pre-storm mean for reference

# Professional color palettes (Economist/NYT style)
period_palette <- c(
  "pre" = "#2C5F7D",      # Muted blue-gray for pre-storm
  "during" = "#C0392B",   # Deep red for during storm
  "post" = "#E67E22"      # Warm orange for post-storm
)

p <- ggplot(panel,
  aes(x = rel_time,
      y = effort_n_vessels_km2)) + 
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "#4A4A4A",
             linewidth = 0.5) +
  geom_hline(yintercept = c(period$mean[period$period == "pre"]),
             color = "#4A4A4A",
             linewidth = 1,
             linetype = "dotted") +
  stat_summary(geom = "pointrange", 
               fun.data = mean_se, 
               aes(color = period),
               linewidth = 0.8,
               size = 0.6) +
  theme_minimal(base_family = "serif") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "#4A4A4A", linewidth = 0.5),
    axis.ticks = element_line(color = "#4A4A4A", linewidth = 0.5),
    axis.text = element_text(size = 9, color = "#4A4A4A"),
    axis.title = element_text(size = 10, color = "#2C2C2C", face = "bold"),
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.95),
    legend.justification.inside = c(0, 1),
    legend.title = element_text(size = 10, face = "bold", color = "#2C2C2C"),
    legend.text = element_text(size = 9, color = "#4A4A4A"),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15),
    text = element_text(color = "#2C2C2C"),
    panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Days relative to first / last day with storm-force winds (> 18 m/s)",
       y = quote("Vessel activity (# vessels/"~km^2~")"),
       color = "Period") +
  scale_color_manual(values = period_palette) +
  scale_x_continuous(breaks = seq(-60, 45, by = 15))


## Export: save the figure to the results directory. File name chosen to
## match other project outputs and be easy to reference in reports.
ggsave(plot = p,
       filename = "results/img/ts_effort.png",
       width = 8,
       height = 4)
