################################################################################
# title
################################################################################
#
# Juan Carlos Villaseñor-Derbez
# jc_villasenor@miami.edu
# date
#
# Description
#
################################################################################

# SET UP #######################################################################

## Load packages ---------------------------------------------------------------
pacman::p_load(
  here,
  tidyverse,
  rnaturalearth,
  sf
)


## Load data -------------------------------------------------------------------
### Hurricane data
#### Column names
col_names <- names(read_csv("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.NA.list.v04r01.csv", n_max = 0))
#### Data
hur <- read_csv("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.NA.list.v04r01.csv",
                col_names = col_names,
                skip = 2) |> 
  janitor::clean_names()

# Gulf states
gulf_states <- ne_states(country = c("United States of America")) |> 
  filter(name %in% c("Florida", "Alabama", "Mississippi", "Louisiana", "Texas"))

# sf_use_s2(F)
# mex_us <- ne_countries(scale = "large") |>
#   st_crop(st_buffer(gulf_states, 2.5))
# sf_use_s2(T)

# Build stats table for the Gulf. These come from the spreadhseet put together by Amelia.
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
                      "Jobs: ", format(employment, big.mark = ",", scientific = FALSE), "\n",
                      "Revenue: ", format(income, big.mark = ",", scientific = FALSE), "\n",
                      "Value: ", format(value_added, big.mark = ",", scientific = FALSE), "\n",
                      "Production: ", format(production, big.mark = ",", scientific = FALSE), "\n"))

# PROCESSING ###################################################################

## Identify relevant hurricanes ------------------------------------------------
gulf_hurs <- hur |> 
  filter(!name == "UNNAMED", # Retain only named storms
         season >= 2000) |>  # Retain only data since 2000
  select(name, year = season, date = iso_time,
         lon, lat,
         wind = usa_wind, sshs = usa_sshs)

# Make a spatial version of hurricanes
hur_sf <- gulf_hurs |> 
  st_as_sf(coords = c("lon", "lat"),
           crs = "EPSG:4326") |> 
  group_by(name, year) |> 
  summarize(max_sshs = as.character(max(sshs)),
            .groups = "drop",
            do_union = F) |> 
  st_cast("LINESTRING") |> 
  st_filter(gulf_states)

# Perform a spatial join
gulf_states_affected <- st_join(gulf_states |>
                                  rename(state_name = name), 
                                hur_sf) |> 
  count(state_name) |> 
  left_join(gulf_stats, by = join_by("state_name" == "state")) |> 
  mutate(state_name = fct_relevel(state_name, "Texas", "Louisiana", "Mississippi", "Alabama", "Florida")) |> 
  arrange(state_name)


# VISUALIZE ####################################################################

## Another step ----------------------------------------------------------------
p1 <- ggplot() + 
  geom_sf(data = gulf_states_affected,
          mapping = aes(fill = n),
          color = "#4A4A4A",
          linewidth = 0.5) +
  geom_label(data = gulf_states_affected, 
             x = c(-107.5, -94.5, -93, -85, -87),
             y = c(34, 27, 35, 34, 27),
             mapping = aes(label = text),
             hjust = "left",
             color = "black",
             fill = "white",
             alpha = 0.5,
             label.r = unit(0.1, "lines"),
             size = 3,
             family = "serif") +
  theme_void(base_family = "serif") +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.11),
        legend.direction = "horizontal",
        legend.background = element_rect(color = "black"),
        legend.margin = margin(t = 15, b = 10, l = 10, r = 10),
        legend.title = element_text(size = 10, face = "bold", color = "#2C2C2C"),
        legend.text = element_text(size = 9, color = "#4A4A4A"),
        plot.margin = margin(10, 10, 10, 10),
        text = element_text(color = "#2C2C2C")) +
  scale_fill_gradient(low = "#E8F4F8",
                      high = "#8B3A3A",
                      limits = c(0, 50),
                      name = "Number of\nHurricanes") +
  labs(x = NULL,
       y = NULL)

p1
# EXPORT #######################################################################

## The final step --------------------------------------------------------------
ggsave(plot = p1,
       filename = "results/img/map.png",
       width = 8,
       height = 4.5)

