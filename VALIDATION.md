# Snapshot verification

## Current supported version

The current requirements pin Generapy to
`bd50a507ff009127768dfdf53e209e6c0a364142`
(`0.1.1.dev35+gbd50a507f`). A clean installation, all 17 check scripts and
both complete numerical examples passed on 16 November 2026. The current
code and requirements are intended to be used together. The full
[upgrade validation](#generapy-upgrade-compatibility-test--16-november-2026)
below records the environment, API changes, exact baseline comparison and
remaining work. Earlier sections describe historical checks, rather than
the current installation status.

## Original baseline

Prepared and checked on 15 November 2026 using the existing pinned
environment: Python 3.14.7, Generapy `0.1.1.dev4+gaed252bd5`, SymPy
`1.15.0.dev0`, mpmath 1.4.1 and NumPy 2.5.3 on macOS ARM64.

This verifies the standalone directory layout and copied code in that
environment. At that point, a fresh dependency installation had not yet
been tested; the Generapy upgrade was tested subsequently as recorded below.

## Completed checks

- All 17 check/diagnostic scripts listed by `code/run_checks.py` completed
  successfully. The 16 original suite entries were run together; the
  additional Mumford consistency check was run separately. Some diagnostic
  scripts report residuals for inspection rather than assert tolerances.
- Fresh Generapy geometry and mode evaluation for both examples at
  `t=0` and `t=0.1` agreed with Hamilton's equations. Maximum scaled
  vector-field residuals were approximately `4.2e-16` (hyperelliptic)
  and `1.3e-14` (trigonal). Fresh initial modes agreed with the stored
  initial data.
- The parameter/validation table generator reproduced the included
  manuscript tables from the stored JSON records.
- The staged LaTeX source built independently to 51 pages, with no
  undefined references, multiply defined labels or overfull boxes.
- The manuscript packager collected 25 required source/bibliography/
  figure files. No historical workspace directory is needed to build it.
- All staged Python sources parse, and checks import no helpers from
  the previous repository layout.

The full 41-point analytic comparison plots and long 35-digit RK4 runs
were not repeated during this copy-only preparation. Their existing
parameter records, trajectories, figures and convergence data are
included unchanged in meaning. Check scripts may regenerate their own
small result JSON files.

## Public repository metadata update

On 15 November 2026, author details and the final public repository URL
were added to the manuscript, and the README recorded the supplied CC0
1.0 Universal dedication. The manuscript was rebuilt with `latexmk` to
51 pages, with no undefined references or overfull boxes. The rendered
title page and repository-link paragraph were visually checked.
No numerical code or stored results were changed, and numerical suites
were not rerun for this metadata-only update. Numerator-zero-safe BA
differentiation, the Generapy release migration, and a permanent manuscript
identifier remain open in the README checklist.

## Generapy upgrade compatibility test — 16 November 2026

The `generapy-020` branch was tested against Generapy commit
`bd50a507ff009127768dfdf53e209e6c0a364142`, verified against GitHub's
`main` on the test date. The installed version is
`0.1.1.dev35+gbd50a507f`. The comparison baseline is this project's commit
`206ad552d927531dcb71efb2010ea3867dbc97d4`, using Generapy
`0.1.1.dev4+gaed252bd5`.

### Environment and compatibility change

A separate `.venv-generapy/` environment was created with Python 3.14.7
and installed from `requirements.txt`. Generapy was built and installed
from its pinned GitHub revision, rather than imported from an editable
checkout. Installed VCS metadata confirms the exact commit, and `pip check`
reports no broken requirements. SymPy, mpmath, NumPy and Matplotlib retain
the baseline versions; the original `.venv/` was left unchanged.

The first full-suite attempt found that `kleinian_baker_akhiezer` had been
replaced by the public `baker_akhiezer` API. The two affected scripts,
`latex_final_ba_conventions.py` and `quartic_hyperelliptic_direction.py`,
now call that API with explicit `normalization='spectral'`. The old and
new formulas use the same spectral normalization. The paper's reference
evaluator and numerical-example algorithms were not changed.

### Completed checks and results

- All 17 check scripts completed: the 13 quick/symbolic entries passed
  during the initial full-suite run; the two sigma checks passed after
  the API update; both diagnostics passed. Diagnostic output was reviewed,
  including the fifth-trace result of about `1.5e-14` after 40-digit
  interpolation and al-family Hamilton residuals below `6e-40`.
- Both fresh `--check-only` examples passed at 24 digits. The reported
  `t=0.1` canonical residuals were `4.1e-16` (hyperelliptic) and `1.2e-14`
  (trigonal).
- Both full examples passed on `0 <= t <= 6`, with 41 direct analytic
  samples, 601 plotted ODE samples, and 6000/12000-step RK4 integrations.
- Every array in both saved NPZ files is exactly equal to the baseline:
  analytic/ODE time grids, BA/ODE trajectories, weights, phases and `z0`.
  Geometry JSON differs only in Generapy version; example-summary JSON
  differs only in version and timings. All numerical metrics are identical.
- Hyperelliptic/trigonal maximum canonical residuals remain `1.2e-15` /
  `1.2e-12`; fine ODE discrepancies remain `1.7e-11` / `1.8e-11`.
- Regenerated PNGs are byte-identical. Regenerated figure PDFs are identical
  after excluding their creation dates; the existing PDFs were retained
  to avoid date-only binary changes.
- The spectral-image check passed again against the regenerated data.
  Regenerated manuscript parameter and validation tables are unchanged.

The long 35-digit convergence integrations were not rerun: every stored
input and reference array they consume is exactly unchanged. Their existing
records remain applicable. No manuscript content changed, so the existing
compiled manuscript was retained. Timings are recorded, but these runs
were not a controlled performance benchmark.

### Reproduction commands

From the repository root, with the new requirements pin:

```sh
.venv/bin/python -m venv .venv-generapy
.venv-generapy/bin/python -m pip install -r requirements.txt
.venv-generapy/bin/python -m pip check
.venv-generapy/bin/python code/ba_numeric_examples.py --case both --check-only
.venv-generapy/bin/python code/run_checks.py --suite all --log-dir .check-logs/generapy-candidate
MPLBACKEND=Agg .venv-generapy/bin/python code/ba_numeric_examples.py --case both
.venv-generapy/bin/python code/prescribed_curve_spectral_image.py
.venv-generapy/bin/python code/generate_ba_validation_table.py
```

Per-script logs and the complete-example log from this run are under
`.check-logs/generapy-candidate/` (ignored by Git). The current tracked
example summaries record the candidate version and timings.

This establishes compatibility of the existing calculations with the
unreleased commit. Migration to directional sigma jets or the public
reference-normalized BA jets, and explicit numerator-zero tests, remain
separate work; the zero-safe differentiation task is still open. Testing
an eventual PyPI release remains a later step. No Generapy release or tag
was created.
