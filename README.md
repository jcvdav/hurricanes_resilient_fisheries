# Hurricanes & Resilient Fisheries

Code and data for **"Competitiveness in American Seafood Requires Promoting Climate-Resilient Fisheries"**.

A single R script (`scripts/01_exposure.R`) produces one combined two-panel figure characterizing hurricane exposure across Gulf-state coastal counties, plus a LaTeX table of Gulf-state fisheries statistics.

---

## Outputs

### Figure — `results/img/annual_and_total_exposure.png`

A two-panel composite (panels labeled A and B).

**Panel A — Annual hurricane-season county-days exposed.** Bars show the annual sum of county-days exposed to storm-force winds (sustained winds ≥ 17.5 m/s) during hurricane season (Jun–Nov), restricted to coastal counties in the five Gulf states (Texas, Louisiana, Mississippi, Alabama, Florida). The overlaid line + points is a 3-year running average computed over consecutive calendar years (zero-filled for years with no qualifying exposure so the window slides over consecutive years).

**Panel B — Gulf Coast hurricane exposure map.** A choropleth of the five Gulf states shaded by the number of distinct (year, hurricane) pairs that produced storm-force sustained winds in any sampled county of that state.

### Table — `results/tab/fisheries_stats.tex`

A LaTeX table of Gulf-state commercial and recreational fisheries statistics (Employment, Production, Income, Sales, Value Added), split by sector. Generated with `modelsummary::datasummary` and post-processed to (1) strip `siunitx` `\num{...}` wrappers so comma-formatted numbers render verbatim, and (2) wrap the table in `\begingroup\footnotesize ... \endgroup` to fit on a portrait page.

---

## Data sources

- **County-day sustained-wind series:** `data/raw/hurricanes/2024_forecast_data.csv` — per hurricane / county / day sustained-wind reconstructions.
- **Sample counties (analysis panel FIPS):** `data/processed/sample_counties.csv` — filtered to the five Gulf states; further intersected with the Gulf of Mexico geometry (MarineRegions `mrgid = 4288`, buffered 50 km) via `mregions2` and `tigris` county shapes.
- **Fisheries statistics:** `data/raw/noaa_stats/gulf_fisheries_stats.csv` — Gulf-state estimates compiled from NOAA's Fisheries Economics of the U.S. report.

---

## Repository structure

```
.
├── Makefile                              # `make` builds the figure
├── scripts/
│   └── 01_exposure.R                     # Produces the figure and the table
├── data/
│   ├── raw/
│   │   ├── hurricanes/2024_forecast_data.csv
│   │   └── noaa_stats/gulf_fisheries_stats.csv
│   └── processed/
│       └── sample_counties.csv
└── results/
    ├── img/annual_and_total_exposure.png
    └── tab/fisheries_stats.tex
```

---

## Usage

```bash
make           # build the figure (and table)
make clean     # remove generated figure(s)
```

Or run the script directly:

```bash
Rscript scripts/01_exposure.R
```

---

## Dependencies

R packages (loaded via `pacman::p_load`): `here`, `tidyverse`, `lubridate`, `rnaturalearth`, `sf`, `tigris`, `cowplot`, `mregions2`, `janitor`, `modelsummary`.

---

## Author

Juan Carlos Villaseñor-Derbez — [jc_villasenor@miami.edu](mailto:jc_villasenor@miami.edu)
