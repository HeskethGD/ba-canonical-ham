# Baker–Akhiezer functions in modified Jiang–Rauch-Wojciechowski systems

Manuscript and reproducibility code for *Baker–Akhiezer functions on
arbitrary curves in modified Jiang–Rauch-Wojciechowski systems*.

The paper studies a canonical cubic multiwave Hamiltonian, its complete
integral count and spectral curves, and the two-way correspondence between
its vector–matrix–vector case and scalar Baker–Akhiezer (BA) data.
The higher spectral Hamiltonians generate Abel translations with explicit
constant mode phases. Numerical examples use genus-two hyperelliptic and
genus-three trigonal curves.

Author: Graham Hesketh, Independent Researcher, UK
([gdh1e10@gmail.com](mailto:gdh1e10@gmail.com)).

Public repository: <https://github.com/HeskethGD/ba-canonical-ham>.

## License

The original manuscript, code, figures and numerical records in this
repository are dedicated to the public domain under **CC0 1.0 Universal**;
see [`LICENSE`](LICENSE). External dependencies retain their own licenses.

## Contents

- [`latex/main.pdf`](latex/main.pdf): compiled manuscript.
- [`latex/main.tex`](latex/main.tex): manuscript source, with all inputs,
  bibliography and figures in the same directory tree.
- [`code/README.md`](code/README.md): check inventory and reproduction commands.
- `code/results/`: recorded parameters, BA/ODE trajectories and check outputs.
- `code/figures/`: black-and-white figures generated from the examples.
- `code/support/`: small local helpers needed by the exact checks.
- [`VALIDATION.md`](VALIDATION.md): current validation results, tested
  dependency revisions and historical baseline verification.

The repository contains the material for this paper. Separate determinant
and Kovalevskaya projects, exploratory notebooks, source-paper collections
and review correspondence are not part of this distribution.

## Python environment

The recorded runs used Python 3.14.7 on macOS ARM64. Install into a fresh
environment from the repository root:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip check
.venv/bin/python code/run_checks.py --suite quick
```

The current requirements preserve the tested source commits of Generapy
and SymPy and the recorded numerical dependency versions. Generapy is
`0.1.1.dev35+gbd50a507f`, pinned to commit
`bd50a507ff009127768dfdf53e209e6c0a364142` from its unreleased `main`.
A clean installation, all 17 check scripts and both complete numerical
examples passed on 16 November 2026. All stored BA/ODE arrays and numerical
validation metrics reproduce the original baseline exactly. A tagged
release has not yet been selected. The plotting and exact-check code does
not require Marimo.

Git is required to install the pinned Generapy and SymPy sources. When
updating an existing checkout, rerun the installation and dependency check:

```sh
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip check
```

The current check scripts use Generapy's public `baker_akhiezer` API and
require the updated dependency pin. With another virtual environment,
replace `.venv/bin/python` in the commands below, or set the Makefile's
`PYTHON` variable, for example `make check-all PYTHON=/path/to/venv/bin/python`.

## Reproduce the paper

Run commands from the repository root. The main example script computes
periods, lifts and exponent constants once, uses Generapy sigma jets at
sparse analytic sample points, and independently integrates Hamilton's
equations with RK4:

```sh
.venv/bin/python code/ba_numeric_examples.py --case hyperelliptic
.venv/bin/python code/ba_numeric_examples.py --case trigonal
.venv/bin/python code/prescribed_curve_spectral_image.py
```

Defaults reproduce the recorded interval `0 <= t <= 6`, precision of
24 decimal digits, 41 analytic sample points and 601 plotted numerical
points. The trigonal sigma evaluation is substantially slower than the
hyperelliptic example. To check construction without regenerating plots:

```sh
.venv/bin/python code/ba_numeric_examples.py --case both --check-only
```

The 35-digit convergence calculations use the **same stored binary64
coefficients and initial modes**. They isolate integration arithmetic;
they do not reconstruct the curve at higher precision:

```sh
.venv/bin/python code/ba_convergence_diagnostic.py --multiprecision --hyperelliptic
.venv/bin/python code/ba_convergence_diagnostic.py --multiprecision
.venv/bin/python code/generate_ba_validation_table.py
```

Those commands update files in `code/results/` and the generated parameter
and validation tables in `latex/`. The existing records are included so
the manuscript and checks can be inspected before rerunning the expensive
calculations. Regenerating examples replaces their records and figures.

For exact and structural checks, see the suites in `code/README.md`.
They supplement the manuscript proofs with explicit examples and
symbolic identities; finite rank witnesses are not all-size proofs.

## Build the manuscript

A TeX installation with `latexmk`, BibTeX and the packages listed in
`latex/main.tex` is required. From the repository root:

```sh
cd latex
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Alternatively use `make paper`. The source includes `main.bbl`; the
bibliography database is `references.bib`. To prepare the submission ZIP
after rebuilding:

```sh
.venv/bin/python code/package_manuscript.py
```

This copies current code figures into `latex/figures/` and packages the
recursive TeX input graph, bibliography and figures. It does not include
the Python code or private working material in the manuscript ZIP.

## Publication checklist

- [x] Test the unreleased Generapy commit in a fresh environment and
  regenerate both complete examples. Results reproduce the baseline;
  see `VALIDATION.md` for the scope and commands.
- [ ] Replace the tested Generapy development commit pin with the chosen
  tagged release when available, then verify a fresh installation of that
  release. Generapy publication is managed separately.
- [ ] Make reference-ratio differentiation zero-safe before widening the
  evaluation domain to numerator-zero crossings. In
  `code/reference_sigma_ba.py`, replace logarithmic numerator derivatives
  with direct quotient/product jets that divide only by the nonzero
  denominator and reference values. Verify ordinary and dual jets through
  order three at numerator zeros, agreement away from zeros, and both
  numerical examples. The current evaluator and tested sample charts are
  retained in the current version.
- [x] Add author information, the chosen repository license (CC0 1.0
  Universal), and the final public repository URL to the project;
  author details and the repository link are also in the manuscript.
- [ ] Add the permanent manuscript identifier (for example, an arXiv ID
  or DOI) once assigned. The GitHub repository URL is recorded above;
  no permanent manuscript identifier has been supplied yet.

## Further development

The current version uses the tested Generapy dependency and public spectral
BA API. Adopting directional sigma jets and the public reference-normalized
BA jets is planned separately, together with the numerator-zero tests above.
Numerical example records identify the tested Generapy development version;
their numerical values are unchanged from the original baseline.
