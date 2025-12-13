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

# Load packages ----------------------------------------------------------------
library(here)
library(gfwr)
library(stormwindmodel)
library(tidyverse)
library(sf)
library(terra)

# Define user functions --------------------------------------------------------
# Function to fix the length of timestamp strings from Hurricane data
fix_length <- function(char, desired_length){
  current_length <- str_length(char)
  
  missing <- desired_length - current_length
  
  zeroes <- character(length = length(missing))
  
  zeroes[which(missing >= 1)] <- 0
  
  char <- paste0(zeroes, char)
  
  return(char)
}

# First, define a function that returns a time-series for each pixel
get_hurricane_exposure <- function(track, grid) {
  # Calculate exposure
  exposure <- stormwindmodel::calc_grid_winds(hurr_track = track,
                                              grid_df = grid)
  # Extract the time series
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
  
  # Save the data
  return(x = ts)
}

# Load data --------------------------------------------------------------------
GoM <- read_sf("/Users/jcvd/Library/CloudStorage/Box-Box/01_project_data/central/World_Seas_IHO_v3") |> 
  filter(NAME == "Gulf of Mexico") |> 
  st_simplify(dTolerance = 10000)

# Hurricanes
col_names <- names(read_csv("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.NA.list.v04r01.csv", n_max = 0))
hur <- read_csv("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.NA.list.v04r01.csv",
                col_names = col_names,
                skip = 2)

## PROCESSING ##################################################################
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
    mutate(lat = lat + 0.05,
           lon = lon + 0.05,
           grid_id = paste(lon, lat))
  
  return(effort)
}

rast <- map_dfr(2015:2020,
                my_get_raster)



## VISUALIZE ###################################################################

# X ----------------------------------------------------------------------------


## EXPORT ######################################################################

# X ----------------------------------------------------------------------------

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

# Calculate exposure -----------------------------------------------------------

ocean_grid <- rast |> 
  select(glon = lon, glat = lat) |> 
  distinct() |> 
  mutate(gridid = paste(glon, glat),
         glandsea = F) |> 
  as.data.frame()

ocean_grid_rast <- rast(ocean_grid |>
                          select(1:2),
                        crs = "EPSG:4326") |> 
  cellSize() |> 
  as.data.frame(xy = T) |> 
  rename(lon = x, lat = y) |> 
  mutate(area = area / 1e6) |> 
  as_tibble()


mirai::daemons(20)
# For ocean pixels
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

period <- group_by(panel,
                   period) |> 
  summarize(mean = mean(effort_n_vessels_km2,
                        na.rm = T),
            se = sd(effort_n_vessels_km2, na.rm = T) / sqrt(n()))

p <- ggplot(panel,
       aes(x = rel_time,
           y = effort_n_vessels_km2)) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(period$mean[period$period == "pre"])) +
  stat_summary(geom = "pointrange", fun.data = mean_se, aes(color = period)) +
  theme_minimal() +
  labs(x = "Days relative to first / last day with storm-force winds (> 18 m/s)",
       y = quote("Vessel activity (# vessels/"~km^2~")"),
       color = "Period") +
  theme(legend.position = "inside",
        legend.position.inside = c(0, 1),
        legend.justification.inside = c(0, 1)) +
  scale_x_continuous(breaks = seq(-56, 45, by = 7))


ggsave(plot = p,
       filename = "results/img/ts_effort.png",
       width = 6,
       height = 3)
