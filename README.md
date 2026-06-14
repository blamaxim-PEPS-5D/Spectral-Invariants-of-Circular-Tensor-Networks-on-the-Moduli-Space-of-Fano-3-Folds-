================================================================================
Spectral_Invariants_Full_Verification.py
================================================================================

Complete verification suite for the manuscript:
"Spectral Invariants of Circular Tensor Networks on the Moduli Space of Fano 3-Folds"

Author: Massimiliano Blandino
ORCID: https://orcid.org/0009-0006-3252-4011
Repository DOI: (to be assigned upon release)

================================================================================
DESCRIPTION
================================================================================

This script performs all numerical and symbolic verifications required for the
manuscript. It generates 18 outputs and 7 figures, covering:

1. Closed-form α⁻¹ from continued fraction [14,1,7,3,1,3]
2. Historical CODATA analysis (2006-2022, all quotients ≤ 45)
3. MPS Spectral Convergence (ln λ_max → A_geo + π)
4. Monte Carlo simulation (⟨S⟩, σ, ergodic convergence)
5. Sensitivity scan for τ (invariance of ⟨S⟩)
6. Commutator norm (theoretical: exactly zero)
7. Sensitivity scan for ε (non-zero coupling stability)
8. PF–DS Numerical Equivalence (Minkowski coefficients match ID-69)
9. PF–DS Symbolic Verification (SymPy notebook integration)
10. Sensitivity scan for D (bond dimension 40-50)
11. Statistical test for D=45 (p-value, false positive rate)
12. MIS bridging (H_seq → H_Fano, peaked at constant 45)
13. Scan of 105 Fano families (only ID-69 matches)
14. Cut norm convergence (graph sequence → graphon on S³)
15. Spectral gap (Λ₁/Λ₀, Gegenbauer integrals)
16. Graphon parameters (γ, κ_W)
17. Commutator norm table (numerical verification)
18. Minimal reproducible script (fast verification)

================================================================================
REQUIREMENTS
================================================================================

Python 3.10 or higher with the following libraries:
  - mpmath (>= 1.3.0)
  - numpy (>= 1.24.0)
  - scipy (>= 1.10.0)
  - matplotlib (>= 3.7.0)
  - sympy (>= 1.12)

Install all dependencies with:
  pip install mpmath numpy scipy matplotlib sympy

================================================================================
EXECUTION
================================================================================

Run the script from the terminal:

  python Spectral_Invariants_Full_Verification.py

The script will:
  - Generate 18 numbered outputs printed to the console
  - Save 7 figures as PNG files in the current directory
  - Display a final summary with verification status

Total execution time is approximately 5 minutes on standard hardware.

================================================================================
OUTPUT FILES
================================================================================

Figures generated:
  - historical_codata_quotients.png
  - mps_convergence_error.png
  - monte_carlo_convergence.png
  - sensitivity_tau.png
  - sensitivity_D.png
  - cut_norm_convergence.png
  - spectral_gap_integrands.png

================================================================================
REPRODUCIBILITY
================================================================================

The script uses:
  - Fixed random seed (implicitly via numpy, reproducible across runs)
  - 50 decimal digit precision (mpmath.dps = 50)
  - Deterministic numerical integration (mpmath.quad)

To verify reproducibility, run the script multiple times. The numerical outputs
should be identical to the values reported in the manuscript.

================================================================================
ADDITIONAL FILES IN THE REPOSITORY
================================================================================

  - notebooks/pf_ds_symbolic.ipynb : SymPy notebook for PF–DS symbolic reduction
  - null_model_test.py             : null model test script for D=45 uniqueness
  - cf_sensitivity.py              : continued fraction sensitivity analysis
  - commutator_grid_check.py       : script to compute commutator norms
  - sdps/ and logs/sdp/            : SDP implementations and solver logs
  - minimal_script.py              : minimal reproducible script
  - README.md                      : execution instructions and environment specs
  - requirements.txt               : Python dependencies

================================================================================
LICENSE
================================================================================

This work is licensed under the Creative Commons Attribution 4.0 International
License (CC BY 4.0). You are free to share and adapt the material for any
purpose, provided appropriate credit is given to the original author.

================================================================================
CONTACT
================================================================================

Massimiliano Blandino
ORCID: https://orcid.org/0009-0006-3252-4011
