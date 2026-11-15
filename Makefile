PYTHON ?= .venv/bin/python

.PHONY: check check-all examples tables paper package

check:
	$(PYTHON) code/run_checks.py --suite quick

check-all:
	$(PYTHON) code/run_checks.py --suite all

examples:
	MPLBACKEND=Agg $(PYTHON) code/ba_numeric_examples.py --case both

tables:
	$(PYTHON) code/generate_ba_validation_table.py

paper: tables
	cd latex && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

package: paper
	$(PYTHON) code/package_manuscript.py
