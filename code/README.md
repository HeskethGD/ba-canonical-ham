# Code and check inventory

Run these scripts from the repository root with its Python environment.
Only checks supporting the manuscript are included. Names prefixed
`latex_final_` are retained for traceability; the manuscript directory is
now `latex/`.

Install the root [`requirements.txt`](../requirements.txt) before running
the checks, including when updating an existing environment. The supported
Generapy version is pinned there; the spectral BA checks use its public
`baker_akhiezer` API. See the root [setup instructions](../README.md#installation)
and [validation record](../VALIDATION.md) for the tested versions and results.
The commands below use `.venv/bin/python`; substitute your environment's
Python if it has another name.

## Check suites

```sh
.venv/bin/python code/run_checks.py --suite quick
.venv/bin/python code/run_checks.py --suite symbolic
.venv/bin/python code/run_checks.py --suite sigma
.venv/bin/python code/run_checks.py --suite diagnostics
.venv/bin/python code/run_checks.py --suite multipoint
.venv/bin/python code/run_checks.py --suite all
```

The default `quick` suite checks physical signs, polynomial higher flows,
the spectral map from stored data, Abel-direction formulas and a repeated-
phase BA sector. `symbolic` includes larger exact checks; `sigma` computes
fresh Generapy BA values, plus independent adapter checks at numerator
zeros and regular points. `diagnostics` includes independent integrations
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
| `mcl_jrw_matrix_embedding.py` | Exact restricted-star canonical dictionary and arbitrary-phase three-block reduction of Menyuk–Chen–Lee equation (2.13), supporting the MCL appendix. |
| `physical_hamiltonian_convention_checks.py` | Physical `-i` convention, Hamiltonian finite differences, Lax equations and complex-conjugate real form. |
| `latex_final_spectral_flows_checks.py` | Polynomial higher-flow time matrices and quartic interaction expansion. |
| `jrw_three_group_integral_count.py` | Distinct-phase spectral-evaluation and spanning-tree independence witnesses, supporting the distinct-phase integral appendix. |
| `jrw_repeated_phase_general.py` | Exact finite-field ranks and evaluated Poisson brackets for spectral plus principal-moment families. |
| `jrw_phase_cluster_insertion.py` | Induction insertion blocks, elimination identities and moment increments, supporting the complete count. |
| `real_forms_regular_locus_checks.py` | Signed real-form identities and local node-smoothing coefficients. |
| `latex_final_symmetric_prym_checks.py` | Symmetric-fold tangency and matrix duality identities, supporting the Prym appendix. |
| `ba_numeric_examples.py` | Both full sigma-versus-ODE validation examples, figures and stored trajectories. |
| `reference_sigma_ba.py` | Ordinary/dual sign and time adapter for Generapy's public reference-normalized BA jets; directional sigma jets supply the separate pair-mode potential. |
| `reference_ba_jet_checks.py` | Independent scalar sigma-ratio differentiation through order three: ordinary/dual normalization, nonfinal coordinates, numerator zeros and nearby points in genera one/two, and singular chart rejection. |
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

The `multipoint` suite contains four additional checks:

| Script | Purpose |
| --- | --- |
| `multipoint_BA_network_check.py` | Exact rational L=M=N=2 BA kernel, interpolation, phase conversion and all 24 rates. |
| `multipoint_BA_elliptic_check.py` | Generapy 0.2.0 theta jets: elliptic two-puncture basis, all rates, left/right spectral eigenvectors and normalized bilinear kernel. |
| `multipoint_higher_flow_check.py` | Finite-part BA times and vertex phases against independent shifted cubic/quartic trace fields. |
| `multipoint_repeated_phase_moments.py` | Exact conservation of all three quadratic moments under the cubic and higher trace flows. |
| `generate_multipoint_validation_table.py` | Generates the additional manuscript residual table and divisor values from the stored higher-flow JSON. |

The last script is a table generator rather than a suite check. The complete
suite now has 22 scripts. The multipuncture wrappers use the existing
public theta jets; no changes to Generapy are required. A native convenience
API assembling general multipoint bases is not assumed.

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

### Quartic-to-cubic comparison

Run `python code/run_checks.py --suite ahh` to check the quartic canonical
parent, the arbitrary-size block/time/phase dictionary, the induced and
frozen matrix brackets, and complete-fiber BA convolution. The first
three scripts use exact SymPy algebra; the last uses an independent
55-digit elliptic theta calculation with a complete regular fiber.
Reports are stored in `code/results/ahh_*.json`.

### Genus-two two-puncture construction

`multipoint_BA_genus2_check.py` constructs ordinary and dual reference
sections algebraically on a smooth genus-two curve. Generapy 0.2.0
provides the periods, Abel maps and complete first-order theta jets.
The fixed algebraic data are initialized once; each evaluation returns
all modes and their analytic time derivatives.

With `--numerics`, it compares both the cubic flow and a fourth-trace
flow with independent RK4 integration, monitors Hamiltonian, power
charges and shifted spectral traces, and writes the single-page
black-and-white figure. Use `generate_genus2_multipoint_table.py` to
regenerate its manuscript table from the stored numerical report.
The `.npz` trajectory record contains both sparse analytic samples
and dense RK4 samples; no new BA evaluator is implemented.

The precision check is reproducible with:

```sh
python code/multipoint_BA_genus2_check.py --dps 40 --output code/results/multipoint_BA_genus2_precision40.json
```

### Reconstruction from canonical initial modes

`ba_genus2_converse_check.py` starts from a specified generic complex
(1,1,2) launch. SymPy computes its characteristic curve and checks its
sextic discriminant is square-free. Adjugate contractions recover the
ordinary and dual BA divisors. Algebraic reference sections recover all
ten initial modes before any theta evaluation; public Generapy periods,
Abel maps and theta jets then give the full time evolution.

```sh
python code/ba_genus2_converse_check.py --numerics
python code/generate_converse_validation_table.py
```

The resulting genus is the generic network genus two. This is a converse
test, unlike examples launched from prescribed BA data. The default
command performs the reconstruction and analytic rate checks; `--numerics`
also compares with independent RK4 integration.

`independent_elliptic_theta_check.py` checks public Generapy theta jets
and ratios against mpmath's Jacobi theta implementation. This backend
is for independent validation only and is not used by the BA plots.
