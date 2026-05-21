# Makefile for hurricanes_resilient_fisheries project
# Produces figures in results/img/

all: results/img/annual_and_total_exposure.png

results/img/annual_and_total_exposure.png: scripts/01_exposure.R | results/img
	Rscript scripts/01_exposure.R

results/img:
	mkdir -p results/img

clean:
	rm -f results/img/annual_and_total_exposure.png results/img/map.png

.PHONY: all clean
