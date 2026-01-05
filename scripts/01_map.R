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
col_names <- names(read_csv("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.NA.list.v04r01.csv", n_max = 0))
hur <- read_csv("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.NA.list.v04r01.csv",
                col_names = col_names,
                skip = 2) |> 
  janitor::clean_names()

gulf_states <- ne_states(country = c("United States of America")) |> 
  filter(name %in% c("Florida", "Alabama", "Mississippi", "Louisiana", "Texas"))

sf_use_s2(F)
mex_us <- ne_countries(scale = "large") |> 
  st_crop(st_buffer(gulf_states, 2.5))
sf_use_s2(T)

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

## Some step -------------------------------------------------------------------
gulf_hurs <- hur |> 
  filter(!name == "UNNAMED",
         season >= 2000) |> 
  select(name, year = season, date = iso_time,
         lon, lat,
         wind = usa_wind, sshs = usa_sshs)

hur_sf <- gulf_hurs |> 
  st_as_sf(coords = c("lon", "lat"),
           crs = "EPSG:4326") |> 
  group_by(name, year) |> 
  summarize(max_sshs = as.character(max(sshs)),
            .groups = "drop",
            do_union = F) |> 
  mutate(max_sshs_bin = case_when(max_sshs == "0" ~ "TS",
                                  max_sshs <= 3 ~ "Minor hurricane",
                                  max_sshs >= 4 ~ "Major hurricane")) |> 
  filter(max_sshs >= 0) |> 
  st_cast("LINESTRING") |> 
  st_filter(gulf_states) |> 
  st_crop(st_buffer(gulf_states, 500000))

gulf_states_affected <- st_join(gulf_states |>
                                  rename(state_name = name), 
                                hur_sf) |> 
  count(state_name) |> 
  left_join(gulf_stats, by = join_by("state_name" == "state")) |> 
  mutate(state_name = fct_relevel(state_name, "Texas", "Louisiana", "Mississippi", "Alabama", "Florida")) |> 
  arrange(state_name)


# VISUALIZE ####################################################################
# Professional color palettes (Economist/NYT style)
# Hurricane categories - refined, publication-quality colors
binned_hurricane_palette <- c(
  "#2C5F7D",  # Tropical Storm - muted blue-gray
  "#E67E22",  # Minor hurricane - warm orange
  "#C0392B"   # Major hurricane - deep red
)


## Another step ----------------------------------------------------------------
p1 <- ggplot() + 
  geom_sf(data = mex_us,
          color = "#D3D3D3",
          fill = "#F5F5F5",
          linewidth = 0.3) +
  geom_sf(data = gulf_states_affected,
          mapping = aes(fill = n),
          color = "#4A4A4A",
          linewidth = 0.5) +
  # geom_sf(data = hur_sf,
  #         mapping = aes(color = max_sshs_bin),
  #         alpha = 0.7,
  #         linewidth = 0.8) +
  geom_label(data = gulf_states_affected, 
             x = c(-100, -91, -92, -85, -82),
             y = c(32.25, 27, 36.75, 36.75, 32.25),
             mapping = aes(label = text),
             color = "black",
             fill = "white",
             label.padding = unit(0.5, "lines"),
             label.size = NA,
             label.r = unit(0.1, "lines"),
             size = 3.2,
             family = "serif",
             nudge_x = c(0, -3, 0, 0, 0),
             nudge_y = c(0, 5, 0, 0, 0)) +
  geom_text(aes(x = -100, y = 38, label = "Unietd States"), size = 5) +
  geom_text(aes(x = -103, y = 25, label = "Mexico"), size = 5) +
  geom_text(aes(x = -91, y = 24, label = "Gulf of Mexico"), size = 5) +
  theme_void(base_family = "serif") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 15, b = 10),
    legend.title = element_text(size = 10, face = "bold", color = "#2C2C2C"),
    legend.text = element_text(size = 9, color = "#4A4A4A"),
    plot.margin = margin(10, 10, 10, 10),
    text = element_text(color = "#2C2C2C")
  ) +
  scale_color_manual(
    values = binned_hurricane_palette,
    name = "Hurricane Category",
    guide = guide_legend(
      override.aes = list(linewidth = 1.2, alpha = 0.8),
      order = 1
    )
  ) +
  scale_fill_gradient(low = "#E8F4F8",
                      high = "#8B3A3A",limits = c(0, 50),
                      name = "Number of\nHurricanes",
                      guide = guide_colorbar(
                        barwidth = 12,
                        barheight = 0.6,
                        title.position = "top",
                        title.hjust = 0.5,
                        order = 2
                      )
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(-109.1, -77.55)) +
  scale_y_continuous(expand = c(0, 0), limits = c(22.1, 39)) +
  labs(x = NULL,
       y = NULL)

p1
# EXPORT #######################################################################

## The final step --------------------------------------------------------------
ggsave(plot = p1,
       filename = "results/img/map.png",
       width = 8,
       height = 6.133)

