################################################################################
# Annual hurricane-season county-days exposed to storm-force winds
################################################################################
#
# ポブラシオナル・クマ
#
# Builds Figure 1: annual county-days exposed to sustained winds >= 17.5 m/s
# during hurricane season (Jun-Nov), for the analysis sample of coastal.
#
################################################################################

# SET UP #######################################################################

## Load packages ---------------------------------------------------------------
pacman::p_load(
  here,
  tidyverse,
  lubridate,
  rnaturalearth,
  sf,
  tigris,
  cowplot,
  mregions2
)

options(dplyr.summarise.inform = FALSE)

## Load data -------------------------------------------------------------------
# Get Gulf shapefile
# gaz_search("Gulf of Mexico") #mrgid = 4288
GoM <- gaz_geometry(x = 4288) |> 
  st_union() |> 
  st_buffer(50000)

### Gulf states (restrict the exposure tally to Gulf-of-Mexico counties)
gulf_states <- c("Texas", "Louisiana", "Mississippi", "Alabama", "Florida")

### Sample counties (analysis panel FIPS, filtered to Gulf states)
sample_counties <- read_csv(
  here("data", "processed", "sample_counties.csv"),
  col_types = cols(fips = col_character())) |>
  filter(state %in% gulf_states)

#Using FIPS code '48' for state 'Texas'
# Using FIPS code '22' for state 'Louisiana'
# Using FIPS code '28' for state 'Mississippi'
# Using FIPS code '01' for state 'Alabama'
# Using FIPS code '12' for state 'Florida'
counties <- counties(state = c(48, 22, 28, 01, 12)) |> 
  mutate(fips = paste0(STATEFP, COUNTYFP)) |> 
  select(fips) |> 
  semi_join(sample_counties, by = "fips") |>
  st_transform(crs = "EPSG:4326") |> 
  st_filter(GoM)

### Hurricane forecast data (per hurricane / county / day sustained winds)
forecast_data <- read_csv(
  here("data", "raw", "hurricanes", "2024_forecast_data.csv"),
  show_col_types = FALSE
)

# PROCESSING ###################################################################

## Build annual exposure series ------------------------------------------------
annual_exposure <- forecast_data |>
  mutate(date = as.Date(date),
         fips = sprintf("%05d", county)) |>
  semi_join(counties, by = "fips") |>
  group_by(hurricane, fips, date) |>
  summarise(wind_obs = mean(vmax_sust_obs, na.rm = TRUE)) |>
  group_by(fips, date) |>
  summarise(wind_obs = max(wind_obs, na.rm = TRUE)) |>
  mutate(year = year(date),
         month = month(date),
         exposed = wind_obs >= 17.5) |>
  filter(month %in% 6:11) |>
  group_by(year) |>
  summarise(county_days_exposed = sum(exposed, na.rm = TRUE)) |>
  # Fill calendar years with no qualifying exposure so the 3-yr window slides
  # over consecutive years (not just over years that appear in the data).
  complete(year = full_seq(year, 1), fill = list(county_days_exposed = 0)) |>
  arrange(year) |>
  mutate(rolling_avg_3y = (county_days_exposed +
                           lag(county_days_exposed) +
                           lag(county_days_exposed, 2)) / 3)

# VISUALIZE ####################################################################

## Bars + 3-yr running average -------------------------------------------------
p1 <- ggplot(annual_exposure,
             mapping = aes(x = year, y = county_days_exposed)) +
  geom_col(mapping = aes(fill = "County-days"),
           width = 0.75) +
  geom_line(mapping = aes(y = rolling_avg_3y,
                          color = "3-yr Running average"),
            linewidth = 1.15,
            na.rm = TRUE) +
  geom_point(mapping = aes(y = rolling_avg_3y,
                           color = "3-yr Running average",
                           fill = "3-yr Running average"),
             size = 2.2,
             shape = 21,
             na.rm = TRUE) +
  scale_x_continuous(name = "Year",
                     breaks = seq(min(annual_exposure$year),
                                  max(annual_exposure$year),
                                  by = 2)) +
  scale_y_continuous(name = "County-days exposed to storm-force winds",
                     labels = scales::comma_format()) +
  scale_fill_manual(values = c("3-yr Running average" = "#0072B2",
                               "County-days" = "#D55E00"),
                    aesthetics = c("color", "fill"),
                    name = NULL) +
  guides(color = "none") +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.95),
    # legend.direction = "horizontal",
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.1, "lines")
  )

## Panel B - Map ---------------------------------------------------------------
# Gulf states
gulf_states <- ne_states(country = c("United States of America")) |> 
  filter(name %in% c("Florida", "Alabama", "Mississippi", "Louisiana", "Texas")) |> 
  select(state = name)

# Build stats table for the Gulf. These come from the spreadhseet put together by Amelia.
# Sources are:
#   - Employment, Income, and Value added from: https://media.fisheries.noaa.gov/2024-07/FEUS-2022-v04-0.pdf Page 9, Table 3.
#   - Production comes from: https://www.fisheries.noaa.gov/foss/, using commercial, year = 2022, region type = NMFS regions, all species and reporting results by year/state.


gulf_stats <- tribble(
  ~ "state", ~"employment", ~"income", ~"value_added", ~"production",
  "Alabama",	6971,	139904,	195934,	13872,
  "Louisiana",	32514,	623462,	837860,	413855,
  "Mississippi",	6954,	126921,	165074,	141018,
  "Texas",	41171,	1259222,	2048829,	27889,
  "Florida",	121710,	4578530,	8206789,	29448
) |> 
  mutate_at(.vars = vars(3:4), \(x) round(x/1e3)) |> 
  mutate(text = paste(state, "\n",
                      "Value: ", format(value_added, big.mark = ",", scientific = FALSE), "\n",
                      "Jobs: ", format(employment, big.mark = ",", scientific = FALSE), "\n",
                      "Production: ", format(production, big.mark = ",", scientific = FALSE), "\n",
                      "Revenue: ", format(income, big.mark = ",", scientific = FALSE), "\n"))

# Calculate exposure by state
exposure_by_state <- forecast_data |>
  filter(vmax_sust_obs >= 17.5) |> 
  mutate(date = as.Date(date),
         fips = sprintf("%05d", county)) |>
  inner_join(sample_counties, by = "fips") |>
  mutate(year = year(date)) |> 
  distinct(year, state, hurricane) |> 
  count(state)

gulf_states_affected <- left_join(gulf_states, exposure_by_state, by = join_by(state)) |> 
  left_join(gulf_stats, by = join_by(state)) |> 
  mutate(state = fct_relevel(state, "Texas", "Louisiana", "Mississippi", "Alabama", "Florida")) |> 
  arrange(state)

p2 <- ggplot() + 
  geom_sf(data = gulf_states_affected,
          mapping = aes(fill = n),
          color = "black",
          linewidth = 0.5) +
  geom_label(data = gulf_states_affected, 
             x = c(-107.5, -94.5, -94.5, -85, -87),
             y = c(34, 27, 35, 34, 27),
             mapping = aes(label = text),
             hjust = "left",
             color = "black",
             fill = "white",
             alpha = 0.75,
             label.r = unit(0.1, "lines"),
             size = 4,
             family = "serif") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.1),
    legend.direction = "horizontal") +
  scale_fill_gradient(low = "#E8F4F8",
                      high = "#8B3A3A",
                      name = "Number of\nhurricanes") +
  labs(x = NULL,
       y = NULL)

p <- plot_grid(p1, p2, ncol = 1, labels = c("A)", "B)"))

# EXPORT #######################################################################

## Save figure -----------------------------------------------------------------
ggsave(plot = p,
       filename = here("results/img/annual_and_total_exposure.png"),
       width = 12,
       height = 12,bg = "white",
       dpi = 300)
