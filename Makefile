# Makefile for hurricanes_resilient_fisheries project
# Produces figures in results/img/

all: results/img/map.png results/img/ts_effort.png

results/img/map.png: scripts/01_map.R | results/img
	Rscript scripts/01_map.R

results/img/ts_effort.png: scripts/02_ts_effort.R | results/img
	Rscript scripts/02_ts_effort.R

results/img:
	mkdir -p results/img

clean:
	rm -f results/img/map.png results/img/ts_effort.png

.PHONY: all clean
