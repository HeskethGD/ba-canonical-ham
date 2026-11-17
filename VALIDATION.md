# Snapshot verification

## Current supported version

The current requirements pin `generapy==0.2.0` from PyPI. Its runtime source
matches the previously tested development revision. A clean installation
was verified for the release; see the latest
[release validation](#generapy-020-pypi-release) below.
All 18 check scripts and both complete numerical examples passed with the
release on 16 November 2026. The [migration validation](#public-ba-jet-migration)
records the BA adapter conventions and numerator-zero checks. The code and
requirements are intended to be used together. Earlier records below
describe historical checks, rather than the current installation status.

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

This established compatibility of the existing calculations with the
unreleased commit. At that stage, migration to the public reference-normalized
BA jets and explicit numerator-zero tests remained separate work. Those were
subsequently completed below. Testing an eventual PyPI release remains a
later step. No Generapy release or tag was created.

## Public BA jet migration

Completed on 16 November 2026 using the same pinned Generapy build and
`.venv-generapy/` environment. The comparison baseline is commit
`e1ea551232264f6db83b97fba3b96a33a4e18f2e` (the dependency compatibility
upgrade merged to `main`). No dependency revision changed in this migration.

### Implementation

`ReferenceSigmaBA` now adapts the paper's ordinary/dual signs and geometric
Abel time to Generapy's public `baker_akhiezer_jet`, with reference
normalization and explicit period/lift data for both curves. All BA values
and derivatives come from that API; local logarithmic differentiation and
reference-normalizer assembly have been removed.

Both ordinary and dual calls use the same moving argument and direction.
Negating the ordinary shift and exponent rate preserves its convention
because sigma parity cancels in the reference ratio. The adapter requires
the half-integral characteristic used by the examples. The supplied rate
describes evolution in the selected coordinate, not arbitrary Abel directions.

The separate `GenerapySigma` adapter uses public directional sigma jets for
the potential and its derivative in the pair-mode formulas. It matches the
BA API's outer guard precision so requests for the common denominator can
reuse Generapy's precision-aware theta cache. No additional result cache is
implemented in this repository.

### Validation

- All 18 scripts in `run_checks.py --suite all` pass, including the new
  `reference_ba_jet_checks.py` in the sigma suite. The additional ordinary
  reference-zero case was also verified in a final direct run of that script.
- The adapter is compared with independent numerical differentiation of
  direct scalar sigma ratios through order three. Tests cover odd/even
  characteristics, a nonfinal coordinate, complex displacements, both jet
  orders, and normalization at the reference. Numerator-zero tests in genera
  one and two cover ordinary and dual modes at the zero and offsets of
  `+/-1e-12`; all derivatives remain finite and agree within `1e-32` at
  40-digit working precision. Singular moving denominators and ordinary/dual
  reference normalizers are rejected.
- Both fresh example checks and the complete 41-point, `0 <= t <= 6`
  examples pass. Maximum canonical residuals remain `1.2e-15` (hyperelliptic)
  and `1.2e-12` (trigonal); fine ODE discrepancies remain `1.7e-11` and
  `1.8e-11`.
- Every trigonal NPZ array is exactly equal to the baseline. Hyperelliptic
  BA and ODE arrays differ by at most `2.3e-29` in absolute value; grids,
  weights, phase coefficients and `z0` are exactly unchanged. Example JSON
  records now name the public BA/directional-sigma backend.
- Both complete examples regenerate identical PNGs and identical figure
  PDFs apart from creation dates. Existing PDFs were retained to avoid
  date-only changes.
- The standard trigonal convergence diagnostic and both 35-digit convergence
  calculations were rerun against the new stored trajectories. Refinement
  ratios are `15.999843` (hyperelliptic) and `15.961593` (trigonal), both
  rounding to `16.0` in the manuscript. The recorded convergence values
  reproduce the previous values; the trigonal JSON additionally records
  the derived successive-difference ratio.
- Regenerated manuscript tables are unchanged. The methods paragraph now
  names the public reference BA and directional sigma jets. The manuscript
  rebuilds to 51 pages without undefined references or overfull boxes; the
  rendered methods page and validation table were visually checked.

Logs and the baseline-array comparison are under `.check-logs/public-ba/`
(ignored by Git). The full examples and convergence commands replace their
tracked result records. Reproduce the migration checks with:

```sh
MPLBACKEND=Agg .venv-generapy/bin/python code/ba_numeric_examples.py --case both
.venv-generapy/bin/python code/run_checks.py --suite all --log-dir .check-logs/public-ba
.venv-generapy/bin/python code/ba_convergence_diagnostic.py
.venv-generapy/bin/python code/ba_convergence_diagnostic.py --multiprecision --hyperelliptic
.venv-generapy/bin/python code/ba_convergence_diagnostic.py --multiprecision
.venv-generapy/bin/python code/generate_ba_validation_table.py
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd latex/main.tex
```

The zero-safe differentiation task is complete for this adapter. Sigma
denominator and reference zeros still require a different regular chart.
Timing records are observations from validation runs, not controlled speed
comparisons. Generapy publication remains separate.

## Generapy 0.2.0 PyPI release

Validated on 16 November 2026 against this project's preceding commit
`66f35564828ed225fbcb86fbf34accae164146ec`, which used the public BA jets
with Generapy `0.1.1.dev35+gbd50a507f`. Requirements now use
`generapy==0.2.0` from PyPI; the SymPy source pin and other direct dependency
versions are unchanged.

### Source and installation verification

The local Generapy `v0.2.0` tag resolves to
`cd98830890aae1ae5d231b13980345733775bd2e`. Compared with the previously
tested `bd50a507ff009127768dfdf53e209e6c0a364142`, the only tracked changes
are release notes and documentation configuration. All 25 Python source
files packaged in the published wheel match the release tag byte for byte,
excluding generated `_version.py` metadata, which now identifies `0.2.0`.

The wheel downloaded from PyPI is `generapy-0.2.0-py3-none-any.whl`, SHA256:

```text
96263c3664fc8a6b770c13a966fb801b05819adeade6e20d29de71816afde810
```

A fresh `.venv-generapy-release/` environment was created with Python
3.14.7 and installed from the updated requirements. Generapy imports from
that environment's site-packages, reports `0.2.0`, and is not a VCS/editable
installation. `pip check` passes. The prior environments were retained.

### Project validation

- All 18 check scripts pass, including ordinary/dual numerator-zero checks.
- Both complete 41-point examples pass over `0 <= t <= 6`.
- Every saved NPZ array is exactly equal to the preceding development-build
  baseline, including BA/ODE trajectories, weights, phases and time grids.
  Geometry JSON differs only in Generapy version; example summaries differ
  only in version and timings. All numerical validation metrics are identical.
- Regenerated PNGs are identical; figure PDFs differ only in creation dates,
  so existing PDFs were retained. Regenerated manuscript tables are unchanged.

The convergence integrations were not repeated for this release-only change:
their stored inputs and references are exactly unchanged, and those calculations
do not call Generapy. The existing manuscript PDF was retained because its
source, figures and tables did not change. The original code-to-release
comparison and installation/check logs are recorded here and under
`.check-logs/generapy-020-release/` (ignored by Git).

To reproduce from the repository root:

```sh
python3 -m venv .venv-generapy-release
.venv-generapy-release/bin/python -m pip install -r requirements.txt
.venv-generapy-release/bin/python -m pip check
.venv-generapy-release/bin/python code/run_checks.py --suite all
MPLBACKEND=Agg .venv-generapy-release/bin/python code/ba_numeric_examples.py --case both
.venv-generapy-release/bin/python code/generate_ba_validation_table.py
```


## Multipuncture generalization (17 November 2026)

The manuscript now treats all positive L,M,N. Its explicit hyperelliptic,
Mumford and constrained-KP specializations remain one-point results.
Four additional scripts form the `multipoint` suite; `all` now runs 22
checks. No Generapy source was modified.

A fresh ignored `.venv-generapy-release/` was installed with the existing
requirements: Generapy 0.2.0, pinned SymPy 1.15.0.dev0, mpmath 1.4.1,
NumPy 2.5.3 and Matplotlib 3.11.2. `pip check` and all 22 suite scripts
passed. The new elliptic wrapper uses only public Generapy theta jets
with fixed data and a precision-fixed complete-jet cache.

At 55 digits, the L=M=N=2 elliptic example verifies all 24 analytic
canonical rates, both spectral eigenvectors, and the reference-normalized
kernel, with maximum scaled residual below 3e-55. The shifted cubic and
quartic trace checks, with spectral shift 0.17+0.09i, have residuals below
7e-55, including nonzero third-group phase corrections. Their finite-part
rates are constructed independently rather than fitted to the fields.
Exact rational and repeated-phase moment checks also pass.

The numerical geometry has two punctures and degree-two ordinary/dual
divisors on an elliptic curve: it is a genus-one special level of a
network whose generic genus is seven. It does not numerically construct
genus-seven theta data. The general proof is in the manuscript.

The new parameter and residual tables are generated from stored JSON by
`code/generate_multipoint_validation_table.py`. Previous long-time
hyperelliptic/trigonal trajectory and convergence records are retained:
their BA evaluator, coefficients and integration algorithms are unchanged.
The full suite includes their stored spectral-image and BA convention
checks; the long trajectory/convergence calculations were not repeated.

Reproduce the new checks from the repository root:

```sh
.venv-generapy-release/bin/python code/run_checks.py --suite multipoint
.venv-generapy-release/bin/python code/generate_multipoint_validation_table.py
```

The complete run logs are in `.check-logs/multipoint-release/` (ignored).

The rebuilt manuscript has 56 pages (54 before the extension), with no
undefined references or overfull boxes. The 33-file source archive also
rebuilds independently to 56 pages. Reviewed monochrome figure PDFs are
retained when source copies differ only in creation/modification dates.

## Adams–Harnad–Hurtubise comparison

The `ahh` suite checks the selected quartic pullback and its bilinear
projection exactly, verifies the cubic Hamiltonian and phase dictionary
for unequal group sizes including (2,2,2), and derives both matrix
Poisson brackets from their stated definitions. An independent 55-digit
elliptic residue check verifies all six ordered off-diagonal identities
and all three diagonal limits on a complete three-point fiber. The
arbitrary-curve proof is the residue calculation in the manuscript; the
elliptic check is a numerical validation of its signs and frames.

## Genus-two multipuncture flows

The L=M=N=2 genus-two example uses distinct finite punctures and
three-point ordinary/dual BA divisors. The 32-digit checks compare all
24 BA rates with the canonical matrix field, both left and right Lax
eigenvectors with the prescribed curve map, and a fourth-trace flow
with its independently assembled matrix field. An additional 40-digit
run records the precision dependence of these analytic residuals.

The two Abel directions have determinant approximately -0.17423333.
Thus this example tests two independent directions, unlike the elliptic
multipuncture example. Its genus-two curve is a special spectral level
of the generic genus-seven network, with defect five.

The cubic trajectory uses 8000/16000 RK4 steps over [0,8]; the
fourth-trace trajectory uses 32000/64000 steps over [0,1.5]. Each is
compared with 65 direct BA samples. The figure shows 801 dense RK4
samples per flow. Full results and step-halving comparisons are in
`code/results/multipoint_BA_genus2_numerics.json`.
