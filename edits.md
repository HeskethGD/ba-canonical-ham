# Notes for editing paper

Initial clarity and attribution pass completed on 22 November 2026. Checked items have been addressed; unchecked items remain deferred. The general KP correspondence and separate one-point subsection were subsequently incorporated on 22 November 2026, with normalization and flow-matching details in an appendix. The fuller elliptic example remains deferred.

## Section 4 - Correspondences and explicit specializations

- [x] The opening para suggests only one-point exampels are calculable but thats not the case, we clauclate other exampels too, in fact we can calculate all exampels.
  - **Status, 22 November 2026:** Reworded the opening to describe particularly explicit specializations, without restricting calculability to them.

- [x] 4.1 subsection header shoudl perhaps refer to the auhtors and call ti a correspondence, might be a long title though. not sure if good idea.
  - **Status, 22 November 2026:** Renamed to "Correspondence with quartic moment-map systems"; authors remain identified in the opening text.

- [x] "At the function level, a weighted sum of products of two intermediate BA evaluations collapses
  to a single endpoint BA function." - where and when, do we mean as the quartic is rewritten as a cuic system, if os just add a couple fo words to that effect.
  - **Status, 22 November 2026:** The sentence now explicitly locates the collapse in the quartic-to-cubic correspondence.

- [x] we need to add the generalization of our KP correspondence to arbitrary L,M,N. We have this done but not added i think.
  - **Status, 22 November 2026:** Incorporated the arbitrary-N matrix constrained-KP operators, source hierarchy and levelwise canonical generators. Longer normalization and phase-matching arguments are in a dedicated appendix. The genus-two second-time check, stored results and generated table are included in the public repository. A KP field-Poisson equivalence is not claimed.

- [x] the one point specialization stuff that currently sits in the KP section is probably still needed but shoudl come out of KP and have its own section after KP but before hyperelliptic polynomial interpolation.
  - **Status, 22 November 2026:** Split and reordered: the general matrix constrained-KP subsection now precedes a separate one-point specialization, followed by hyperelliptic polynomial interpolation. Preserved the scalar conserved-coefficient calculation with the one-point formulas.

- [x] one points specialization has i clashes "Elementary modes and their canonical equations"
  - **Status, 22 November 2026:** Replaced the first-group index i by ell throughout the one-point formulas, retaining i only as the imaginary unit.

- [x] in "Curve and function conventions" we talk about generapy's conventions as being definitive, but we shoudl not, we shoudl attribute them to source and we can also mention that generapy also uses them. do not refer to specific function names though, those names are now out of date.
  - **Status, 22 November 2026:** Attributed the Baker basis and sigma/BA conventions to mathematical sources, retained agreement with Generapy, and removed API function names.

- [x] "Hyperelliptic BA modes and polynomial interpolation" is it only N=1 that we knwo about polynomial interpolation or shoudl we be discussing this for arbitrary L,M,N?
  - **Status, 22 November 2026:** Clarified the current scope: the compact polynomial formulas use N=1 and the odd-degree infinity marking; the general meromorphic criterion still covers arbitrary L,M,N. A general multipuncture polynomial formula is not claimed.

- [x] citation needed on Weierstrass’s hyperelliptic al function (Weierstrass and Matsutani)
  - **Status, 22 November 2026:** Added Weierstrass (1856) and Matsutani--Previato at the al-function identification, with the original publication and a digitised offprint in the bibliography.

- [ ] "Genus one and many elementary modes" - i think there is an opportunity here ot write a general elliptic funcitons system stripped back nicely to show the pure elegance of the idea on display, use the y^2=4x^3 - g2x -g3 curve deifn points P_j(mu_j) = (\wp(mu_j), \wp(mu_j)') then use mu_j in teh zeta and sigmas etc. Show how this generates a coupled system for arbitrary numbers of genus 1 BA, do it in such a way that it does not have additional nomralization constants floating about, give the value of H in temrs of wp etc make it look liek a pure elliptic funciton identity. What happens in the limit that number of modes tends to infinity? is it a well known lattice model or somehign?
  - **Status, 22 November 2026:** Pending: develop the fuller finite elliptic example separately. The infinite-mode limit also requires a limiting/scaling prescription.

## Section 5 - Numerical Validation

- [x] why do we talk about Green stuff in Numerical Validation, what is that about and what is it used for?
  - **Status, 22 November 2026:** Explained that the bilinear formula evaluates trigonal pair-composite modes and their analytic derivatives independently of Hamilton's equations. Retained the formula and named Green's identity in its derivation.

- [x] "The trigonal example uses genuinely complex weights and phases" put symbols next to these terms to help identify whihc is which.
  - **Status, 22 November 2026:** Added kappa_k and b_k beside weights and phase coefficients.

- [x] "As a separate implementation check, Generapy’s elliptic theta jets through order three,
  theta ratios and logarithmic rates were compared with mpmath.jtheta, including period
  shifts and a point near a theta zero. The maximum scaled discrepancy was below 10−54 at
  55 digits. This independent check is in code/independent_elliptic_theta_check.py; the
  manuscript’s BA evaluations and plots continue to use Generapy." this is silly, generapy has been compared extensively to other systems so this sounds lame, also i work on mpmath anyway so its not really independent. suggest removing and just refer to generapy's info on its validation.
  - **Status, 22 November 2026:** Removed this paragraph and replaced it with a reference to Generapy's documented identity tests and published-result reproductions. The checking script remains in the repository.

## Appendix A

- [x] "The broader wave framework therefore already permits these phase-dependent equations.
  The independent phase parameters are a constant-diagonal extension relative to the printed
  1991 Hamiltonian, as detailed in Appendix B. Their determination from prescribed BA data is
  the role they play in the present paper. The explicit complete-integral calculation in the MCL
  paper concerns the restricted shared-mode family; the all-size canonical count used here is
  stated in theorem 2.2 and proved in Appendix D." - not quite getting all this "The broader wave framework therefore already permits" maybe rephrase this to remove already as no context for what is meant. Are we saying what is differnet about what we do vs what thye do? Lets disucss that.
  - **Status, 22 November 2026:** Rephrased the explicit embedding into the Menyuk--Chen--Lee matrix-wave framework; kept the distinction between their restricted complete-integral calculation and our arbitrary-size count.

# Refs

- [x] Replace arxivs with peer reviewed publications where possible.
  - **Status, 22 November 2026:** First-pass publication audit completed. Most entries already cite published journals/books and retain open author links. No published replacement was verified for Kostov1999, HeskethFWM2026 or HeskethHyperelliptic2026, so those preprints remain. Added the published Weierstrass reference.
