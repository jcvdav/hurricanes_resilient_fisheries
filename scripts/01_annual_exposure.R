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
  lubridate
)

options(dplyr.summarise.inform = FALSE)

## Load data -------------------------------------------------------------------
### Gulf states (restrict the exposure tally to Gulf-of-Mexico counties)
gulf_states <- c("Texas", "Louisiana", "Mississippi", "Alabama", "Florida")

### Sample counties (analysis panel FIPS, filtered to Gulf states)
sample_counties <- read_csv(
  here("data", "processed", "sample_counties.csv"),
  col_types = cols(fips = col_character())
) |>
  filter(state %in% gulf_states)

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
  semi_join(sample_counties, by = "fips") |>
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
  mutate(rolling_avg_3y = zoo::rollmean(county_days_exposed,
                                        k = 3,
                                        fill = NA,
                                        align = "right"))

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
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.1, "lines")
  )

p1

# EXPORT #######################################################################

## Save figure -----------------------------------------------------------------
ggsave(plot = p1,
       filename = "results/img/annual_hurricane_season_county_days_exposed.png",
       width = 12,
       height = 6,
       dpi = 300)
