# Code and check inventory

Run these scripts from the repository root with its Python environment.
Only checks supporting the manuscript are included. Names prefixed
`latex_final_` are retained for traceability; the manuscript directory is
now `latex/`.

## Check suites

```sh
.venv/bin/python code/run_checks.py --suite quick
.venv/bin/python code/run_checks.py --suite symbolic
.venv/bin/python code/run_checks.py --suite sigma
.venv/bin/python code/run_checks.py --suite diagnostics
.venv/bin/python code/run_checks.py --suite all
```

The default `quick` suite checks physical signs, polynomial higher flows,
the spectral map from stored data, Abel-direction formulas and a repeated-
phase BA sector. `symbolic` includes larger exact checks; `sigma` computes
fresh Generapy BA values. `diagnostics` includes independent integrations
and the al-family example. These extended suites can take longer.
`all` does not regenerate the full 41-point figures or the 35-digit RK4
tables; their explicit commands are in the root README.

The runner stops on a nonzero exit status and keeps per-script output in
`.check-logs/`. Some diagnostic scripts print numerical residuals rather
than assert a tolerance; inspect their output when changing algorithms or
dependencies. The mathematical proofs are in the manuscript.

## Manuscript correspondence

| Script | Purpose |
| --- | --- |
| `latex_final_section2_checks.py` | Symbolic Hamilton equations, matrix bracket, phase gauges, trace identities and the Jiang–Rauch-Wojciechowski dictionary. Uses geometric unit brackets; the physical conversion is checked separately. |
| `physical_hamiltonian_convention_checks.py` | Physical `-i` convention, Hamiltonian finite differences, Lax equations and complex-conjugate real form. |
| `latex_final_spectral_flows_checks.py` | Polynomial higher-flow time matrices and quartic interaction expansion. |
| `jrw_three_group_integral_count.py` | Distinct-phase spectral-evaluation and spanning-tree independence witnesses, supporting Appendix B. |
| `jrw_repeated_phase_general.py` | Exact finite-field ranks and evaluated Poisson brackets for spectral plus principal-moment families. |
| `jrw_phase_cluster_insertion.py` | Induction insertion blocks, elimination identities and moment increments, supporting the complete count. |
| `real_forms_regular_locus_checks.py` | Signed real-form identities and local node-smoothing coefficients. |
| `latex_final_symmetric_prym_checks.py` | Symmetric-fold tangency and matrix duality identities, supporting the Prym appendix. |
| `ba_numeric_examples.py` | Both full sigma-versus-ODE validation examples, figures and stored trajectories. |
| `reference_sigma_ba.py` | Reference-normalized sigma BA values and derivatives shared by those examples. |
| `prescribed_curve_spectral_image.py` | Spectral map checked at 36 point/state combinations per stored example. |
| `ba_convergence_diagnostic.py` | Stored-data RK4 precision/refinement checks; multiprecision mode generates the paper's refinement rows. |
| `generate_ba_validation_table.py` | Generate the parameter and validation tables from recorded JSON. |
| `latex_final_ba_conventions.py` | Genera one through three: normalized hyperelliptic BA product, logarithmic rate and Schrödinger identity. |
| `spectral_abel_direction_check.py` | Shifted higher-flow identities and exact Abel-direction spanning example. |
| `quartic_hyperelliptic_direction.py` | Fourth-flow translation and phase matching on fixed genus-two BA data, including a direct Generapy check. |
| `fourth_trace_phase_normalization.py` | Meromorphic finite-part phases, symbolic genus-one identity and exact all-mode genus-two certificates. |
| `fourth_trace_normalization_audit.py` | Fully repeated-phase branch-point realization and finite-part identity. |
| `higher_flow_evolution_diagnostic.py` | Independent genus-one/two integrations versus BA time and phase predictions, including 40-digit coefficient interpolation. |
| `al_function_checks.py` | Al-family mode equations in genera one through three and discriminant diagnostics in genera one and two. |
| `mumford_core_comparison.py` | Classical Mumford brackets, spectral coefficients and divisor flow, supporting Section 6. This does not assert a Poisson equivalence with the full canonical network. |
| `package_manuscript.py` | Synchronize figures and build the source submission ZIP. |

`support/` contains only the symbolic pencil constructor, the exact
genus-two mode embedding and the numerical matrix reference required
by these checks. It does not import another repository or historical
experiment directory.

## Stored files

- `*_sigma_ba_geometry_24.json`: periods, characteristic, fixed-point lifts
  and finite-part rates.
- `*_sigma_ba.json`: version metadata, parameters and validation summaries.
- `*_sigma_ba.npz`: analytic/ODE samples, initial modes and phase constants.
- `*_convergence_diagnostic.json`: convergence diagnostics used by the table.
- Other JSON files: outputs of the correspondingly named check scripts.

The global sigma scale cancels in the reference-normalized ratios.
The exponent rate, Abel basis and lift still need to be mutually
compatible. Sigma calculations use multiprecision arithmetic, while
the assembled example modes and standard RK4 solver use complex128.
