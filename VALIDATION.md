# Snapshot verification

Prepared and checked on 15 November 2026 using the existing pinned
environment: Python 3.14.7, Generapy `0.1.1.dev4+gaed252bd5`, SymPy
`1.15.0.dev0`, mpmath 1.4.1 and NumPy 2.5.3 on macOS ARM64.

This verifies the standalone directory layout and copied code in that
environment. A fresh dependency installation has not yet been tested;
the planned Generapy release migration remains a separate step.

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
