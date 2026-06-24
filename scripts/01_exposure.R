################################################################################
# Annual hurricane-season county-days exposed to storm-force winds
################################################################################
#
# ポブラシオナル・クマ
#
# Builds Figure 1
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
  mregions2,
  janitor,
  modelsummary
)

options(dplyr.summarise.inform = FALSE)

## Load data -------------------------------------------------------------------
# Fisheries production statistics from NOAA
noaa_stats <- read_csv(here("data/raw/noaa_stats/gulf_fisheries_stats.csv")) |> 
  clean_names() |> 
  select(state:units) |> 
  drop_na(state)

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
           color = "black",
           linewidth = 0.5) +
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
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.95),
    legend.key.size = unit(1.1, "lines")
  )

## Panel B - Map ---------------------------------------------------------------
# Gulf states
gulf_states <- ne_states(country = c("United States of America")) |> 
  filter(name %in% c("Florida", "Alabama", "Mississippi", "Louisiana", "Texas")) |> 
  select(state = name)

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
  mutate(state = fct_relevel(state, "Texas", "Louisiana", "Mississippi", "Alabama", "Florida")) |> 
  arrange(state)

sf_use_s2(F)

p2 <- ggplot() + 
  geom_sf(data = ne_countries(country = c("United States of America", "Mexico"), scale = "large") |> 
            st_crop(st_buffer(gulf_states_affected, 0.5)),
          fill = "gray75",
          color = "gray75") +
  geom_sf(data = gulf_states_affected,
          mapping = aes(fill = n),
          color = "black",
          linewidth = 0.5) +
  geom_sf_text(data = gulf_states_affected,
               mapping = aes(label = n),
               nudge_y = 0.5,
               nudge_x = -0.5,
               size = 6) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.9),
    legend.direction = "horizontal") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradient(low = "#FDF1E6",
                      high = "#D55E00",
                      name = "Number of\nhurricanes",
                      guide = guide_colorbar(ticks.colour = "black",
                                             frame.colour = "black")) +
  labs(x = NULL,
       y = NULL)

p <- plot_grid(p1, p2, ncol = 1, labels = c("A)", "B)"))

## Table of Fisheries statistics
data <- noaa_stats |>
  mutate(value = if_else(str_detect(units, "Thousand"), round(value / 1000), value),
         metric = fct_relevel(metric,
                              "Employment",
                              "Production",
                              "Income",
                              "Sales",
                              "Value Added"
                              ),
         state = fct_reorder(state, value, min, .na_rm = T, .desc = T),
         sector = ifelse(str_detect(sector, "Commercial"), "Com.", "Rec.")) |>
  select(State = state, sector, metric, value)

# Keep `value` numeric — turning it into a character makes datasummary treat it
# as another categorical column, exploding the table to hundreds of columns.
# Format with commas via `fmt` instead.
tex_path <- here("results/tab/fisheries_stats.tex")

datasummary(State ~ unique * value * metric * sector,
            data = data,
            fmt = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE),
            output = tex_path,
            title = "\\label{tab:fishery_stats}Commercial and Recreational fisheries in Gulf States are an important source of livelihoods and food.",
            notes = paste(
              "Data come from \\citep{noaa_feus_2022}.",
              "Each economic metric is split by fishing sector: Com. = commercial fisheries; Rec. = recreational fisheries.",
              "Units: Employment is full- and part-time jobs; Production is in metric tons; Income, Sales and Value Added are in millions of U.S. dollars (M USD).",
              "Empty cells indicate that the metric is not separately reported for that state-sector combination."
            ),
            escape = FALSE)

# tinytable wraps numeric-looking cells in \num{...} for siunitx, which then
# parses the embedded comma as a decimal separator and errors out. Drop the
# wrappers so the comma-formatted strings render verbatim. Also wrap the
# whole table in `\begingroup\small ... \endgroup` to fit on a portrait page.
tex <- readLines(tex_path) |>
  str_replace_all("\\\\num\\{([^}]*)\\}", "\\1")
writeLines(c("\\begingroup\\footnotesize", tex, "\\endgroup"), tex_path)

# EXPORT #######################################################################

## Save figure -----------------------------------------------------------------
ggsave(plot = p,
       filename = here("results/img/annual_and_total_exposure.png"),
       width = 9,
       height = 9,
       bg = "white",
       dpi = 300)

