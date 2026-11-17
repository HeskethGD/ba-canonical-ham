PYTHON ?= .venv/bin/python

.PHONY: check check-all examples examples-multipoint tables paper package

check:
	$(PYTHON) code/run_checks.py --suite quick

check-all:
	$(PYTHON) code/run_checks.py --suite all

examples:
	MPLBACKEND=Agg $(PYTHON) code/ba_numeric_examples.py --case both

examples-multipoint:
	MPLBACKEND=Agg $(PYTHON) code/multipoint_BA_genus2_check.py --numerics

tables:
	$(PYTHON) code/generate_ba_validation_table.py
	$(PYTHON) code/generate_multipoint_validation_table.py
	$(PYTHON) code/generate_genus2_multipoint_table.py
	$(PYTHON) code/generate_converse_validation_table.py

paper: tables
	cd latex && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

package: paper
	$(PYTHON) code/package_manuscript.py
