# Canonical multiwave Baker–Akhiezer systems

Manuscript and reproducibility code for *A canonical multiwave Hamiltonian
realization of Baker–Akhiezer functions on arbitrary curves*.

The paper studies a canonical cubic multiwave Hamiltonian, its integrals
and spectral curves, and its correspondence with scalar Baker–Akhiezer
(BA) data for arbitrary positive array sizes L, M and N. The N=1 member
uses a one-point BA pair; N>1 uses scalar multipuncture BA bases on one
curve. Examples include genus-two hyperelliptic, genus-three trigonal,
and elliptic and genus-two L=M=N=2 realizations. The genus-two
multipuncture example compares two independent commuting BA flows.

**Graham Hesketh** · Independent Researcher, UK ·
[gdh1e10@gmail.com](mailto:gdh1e10@gmail.com)

## Contents

- [Paper](latex/main.pdf) and [LaTeX source](latex/main.tex).
- [Code and check inventory](code/README.md).
- [Recorded calculations](code/results/) and [figures](code/figures/).
- [Validation details](VALIDATION.md).

## Installation

Calculations used Python 3.14.7 on macOS ARM64 with Generapy 0.2.0,
SymPy 1.15.0.dev0, mpmath 1.4.1, NumPy 2.5.3 and Matplotlib 3.11.2.
Exact dependencies are pinned in [requirements.txt](requirements.txt).
Git is required to install the pinned SymPy source revision.

From the repository root:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip check
```

Run all commands below from the repository root.

## Numerical examples

Regenerate the BA-versus-ODE comparisons and check their spectral images:

```sh
.venv/bin/python code/ba_numeric_examples.py --case hyperelliptic
.venv/bin/python code/ba_numeric_examples.py --case trigonal
.venv/bin/python code/prescribed_curve_spectral_image.py
```

The examples compare Generapy's reference-normalized BA jets with independent
RK4 integration of Hamilton's equations. Defaults use `0 <= t <= 6`,
24 decimal digits, 41 analytic sample points and 601 plotted ODE points.
Use `--case both` to run both examples in one command. To check their
construction without regenerating trajectories and plots:

```sh
.venv/bin/python code/ba_numeric_examples.py --case both --check-only
```

Recorded results and figures are included. Regeneration replaces files in
`code/results/` and `code/figures/`; `--check-only` also updates geometry records.

## Checks

Run the quick checks or the complete suite of 22 scripts:

```sh
.venv/bin/python code/run_checks.py --suite quick
.venv/bin/python code/run_checks.py --suite all
```

Individual suites and script descriptions are listed in [code/README.md](code/README.md).
Logs are written to `.check-logs/`. The full suite does not regenerate the
complete example plots or the convergence tables.

## Multipuncture BA checks

The public Generapy 0.2.0 theta jets evaluate the new elliptic basis,
kernel and derivatives. The scripts assemble the fixed divisor and
puncture normalization; they do not add a theta-series implementation.
Sigma jets can evaluate the equivalent ratios when those data are supplied.

```sh
.venv/bin/python code/run_checks.py --suite multipoint
.venv/bin/python code/generate_multipoint_validation_table.py
```

This suite checks all 24 canonical rates for L=M=N=2, the spectral
map and normalized eigenvector kernel, higher-flow time/phase matching,
and exact repeated-phase block moments. Its elliptic curve has genus one;
the generic spectral curve of that network has genus seven.

## Convergence calculations and tables

Run the 35-digit RK4 refinement calculations, then generate the manuscript tables:

```sh
.venv/bin/python code/ba_convergence_diagnostic.py --multiprecision --hyperelliptic
.venv/bin/python code/ba_convergence_diagnostic.py --multiprecision
.venv/bin/python code/generate_ba_validation_table.py
.venv/bin/python code/generate_multipoint_validation_table.py
```

These calculations promote the stored binary64 coefficients and initial modes
to higher integration precision; they do not reconstruct the curves at higher
precision. Results are saved in `code/results/`, and tables are written to `latex/`.

## Build the paper

A TeX installation with `latexmk`, BibTeX and the packages listed in
[main.tex](latex/main.tex) is required.

```sh
make paper
```

This regenerates the manuscript tables and builds `latex/main.pdf`.
To compile only the LaTeX source:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd latex/main.tex
```

To build a self-contained manuscript source archive, including its bibliography
and figures:

```sh
make package
```

The archive is written to `latex/submission-source.zip`. Make targets use
`.venv/bin/python` by default; override it with `PYTHON=/path/to/venv/bin/python`.

## License

The original manuscript, code, figures and numerical records are dedicated
to the public domain under [CC0 1.0 Universal](LICENSE). External dependencies
retain their own licenses.

### Genus-two multipuncture example

```sh
.venv/bin/python code/multipoint_BA_genus2_check.py
MPLBACKEND=Agg .venv/bin/python code/multipoint_BA_genus2_check.py --numerics
.venv/bin/python code/generate_genus2_multipoint_table.py
```

The first command checks all 24 mode rates, both spectral eigenvectors
and an independent fourth-trace Abel direction. The numerical run
compares sparse direct BA evaluations with dense RK4 trajectories for
both flows and writes the black-and-white figure and result records.

### Matrix constrained-KP correspondence

```sh
.venv/bin/python code/matrix_kp_hierarchy_check.py
.venv/bin/python code/generate_matrix_kp_validation_table.py
```

This genus-two, two-puncture check evaluates the second normalized KP
time from BA data and compares the ordinary/dual source and kernel
equations with complete Generapy theta jets. It also matches all 24
canonical mode rates to one fixed combination of cubic, fourth-trace
and gauge flows, and checks the ambient source brackets symbolically.
The finite KP time normalization is levelwise; a KP field-Poisson
reduction is not assumed.

### Converse reconstruction

```sh
.venv/bin/python code/ba_genus2_converse_check.py --numerics
.venv/bin/python code/generate_converse_validation_table.py
```

This genus-two (1,1,2) example starts from independent canonical modes,
reconstructs the spectral curve and both BA divisors, and compares the
recovered BA solution with RK4. Its curve has the generic spectral genus.
