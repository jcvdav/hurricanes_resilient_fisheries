# Makefile for hurricanes_resilient_fisheries project
# Produces figures in results/img/

all: results/img/annual_hurricane_season_county_days_exposed.png results/img/map.png

results/img/annual_hurricane_season_county_days_exposed.png: scripts/01_annual_exposure.R | results/img
	Rscript scripts/01_annual_exposure.R

results/img/map.png: scripts/02_map.R | results/img
	Rscript scripts/02_map.R

results/img:
	mkdir -p results/img

clean:
	rm -f results/img/annual_hurricane_season_county_days_exposed.png results/img/map.png

.PHONY: all clean
