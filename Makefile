# make setup                     once: Python environment (.venv) from requirements.txt
# make all                       download, grades, figures, notebook
# make <step> TRAILS="tmb jmt"   one step, some trails (default: every trail in trails.yaml)
PYTHON ?= .venv/bin/python
TRAILS ?=

.PHONY: setup all download grades figures notebook

all: download grades figures notebook

setup:
	python3 -m venv .venv
	.venv/bin/python -m pip install --upgrade pip
	.venv/bin/python -m pip install -r requirements.txt

download:          # OSM tracks and the FABDEM tiles they cross (skips what's there)
	$(PYTHON) scripts/get_data.py $(TRAILS)

grades:            # data/grades/<trail>.csv.gz and summary.csv
	$(PYTHON) scripts/build_grades.py $(TRAILS)

figures:           # figures/density.png, heatmaps.png, zone-bars.png
	$(PYTHON) scripts/plots.py

notebook:          # notebook/walkthrough.html (needs Quarto and the Tour du Mont Blanc's tiles)
	QUARTO_PYTHON=$(PYTHON) quarto render notebook/walkthrough.qmd
