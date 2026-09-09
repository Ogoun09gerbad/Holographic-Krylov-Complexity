# Reproducibility mapping

Precise mapping: **report section → code → input → output → figure/table in the report**. Every entry below was re-run from a clean copy of the repository on 2026-09-09 under Python 3.14.7 (numpy 2.5.2, scipy 1.18.0, sympy 1.14.0, matplotlib 3.11.1); see `REPOSITORY_AUDIT.md` for the full verification log and honest PASS/WARNING status of each.

---

### Sec. 2, "Universal short-time behaviour" — toy 3×3 example (Fig. `fig:toy-ct`)
- **Code:** none found in the repository.
- **Input:** the $3\times3$ Hamiltonian and $b_1=1$ result are stated analytically in the text; no script reproduces `figures/ct_example.png`.
- **Output:** `figures/ct_example.png` (copied as-is from the original project).
- **Report figure:** Fig. `fig:toy-ct`.
- **Status: NOT VERIFIED.** No generating script survives in the project directory. See `SCIENTIFIC_ISSUES.md` / `REPOSITORY_AUDIT.md`.

### Sec. 4, "Global AdS3: analytic benchmark" — numerical validation
- **Code:** `ads3/ads3_validation.py`
- **Input:** $m=1$, $r_0=1.5$, $J\in\{0,1,2,5\}$ for the figure; five $(m,J,r_0)$ triples for the slope-vs-analytic table (hardcoded in `test_cases`, no CLI arguments).
- **Output:** `ads3/holographic_krylov_validation.png`; stdout table of `alpha_analytic` vs `alpha_numeric` and relative error.
- **Report figure/claim:** Fig. `fig:ads3-validation`; text "the relative deviation was at most $1.0\times10^{-7}$".
- **Status: PASS.** Worst relative error reproduced: $1.04\times10^{-7}$.

### Sec. 5/6, "Top-down six-dimensional SCFT backgrounds" / "Two 6d SCFT quiver geometries"
- **Code:** `quiver_comparison/geometry.py`, `quiver_comparison/routhian.py`, `quiver_comparison/run_quiver_comparison.py`
- **Input:** $N=1$, $P=100$, $r_{\rm UV}=\ln(0.01)$, $\eta_0=50$, $J\in\{0,1,3\}$ (hardcoded parameters, documented with their provenance in the script header and in `quiver_comparison/README.md`).
- **Output:** `quiver_comparison/figures/01_geometry_f1_comparison.png` through `08_parameter_scan.png`; `quiver_comparison/results_table.csv`, `quiver_comparison/parameter_scan.csv`.
- **Report figure/table:** Fig. `fig:quiver-f1-comparison` (uses `01_geometry_f1_comparison.png`); Table `tab:quiver-representative` (naive intercepts $0.0877,0.2630,0.2307,0.6920$).
- **Status: PASS.** All four table values reproduced exactly (to the quoted 4 significant figures).

### Sec. 6, "Numerical comparison" — intercept isolation
- **Code:** `quiver_comparison/run_intercept_plot.py`
- **Input:** reuses the same geometry/routhian parameters as above; no new parameters.
- **Output:** `quiver_comparison/figures/09_naive_minus_routhian_intercept.png`.
- **Report figure:** Fig. `fig:quiver-intercept`.
- **Status: PASS.**

### Sec. 6, "Parameter scan" — robustness across quiver size $P$
- **Code:** `quiver_comparison/run_parameter_family_scan.py`
- **Input:** $P\in\{5,20,50,100,200\}$, $N=1$, $J=3$ (for the scan figure).
- **Output:** `quiver_comparison/figures/09_f1_shape_vs_P.png`, `quiver_comparison/figures/10_scan_vs_P.png`, `quiver_comparison/parameter_family_scan.csv`.
- **Report figure:** Fig. `fig:quiver-P-scan`.
- **Status: PASS.**

### Sec. 6, exact $N$-scaling identity $f_i(\eta;N)=\sqrt N f_i(\eta;1)$
- **Code:** self-check inside `quiver_comparison/geometry.py`, invoked from `run_quiver_comparison.py` stdout ("Geometry self-validation" block).
- **Status: PASS** (verified to machine precision in the run log; the report also gives an independent analytic proof in the text, eq. `N-scaling`).

### Sec. 6, symbolic cross-check of Quiver 2's rank function
- **Code:** `quiver_comparison/sympy_verify_quiver2.py`
- **Input:** none (pure symbolic derivation from the ODE $\alpha''=-81\pi^2 R_2(\eta)$ and boundary/continuity conditions).
- **Output:** stdout only — "PASSED: geometry.py's alpha2(eta) matches an independent SymPy re-derivation exactly."
- **Report claim:** "The corresponding $\alpha_2(\eta)$... was independently reconstructed from the rank function and checked symbolically during the project" (Sec. `sec:quivers`, "Quiver 2").
- **Status: PASS.**

### Sec. 7, "Quadratic black-string-chain background"
- **Code:** `blackstring_quiver/geometry.py`, `blackstring_quiver/run_blackstring_quiver.py` (imports `quiver_comparison/routhian.py` by relative path — the two directories must stay siblings).
- **Input:** reference benchmark $P=6$, $\gamma=32$, $(\beta_0,\ldots,\beta_5)=(22,21,19,17,15,15)$; $J\in\{0,1,3\}$ (Table A), $J\in\{0,1,3,8,15\}$ (Table B/Fig. `fig:bs-J-independence`), $r_{\rm UV}\in\{-4.605,-9.200,-20.000\}$ (Table C).
- **Output:** `blackstring_quiver/figures/01_h4_profile.png`, `02_naive_vs_routhian.png`, `03_J_independent_equilibrium.png`; `blackstring_quiver/tableA_mechanism.csv` through `tableD_minimal_quiver.csv`.
- **Report figure/claim:** Fig. `fig:blackstring-h4`, Fig. `fig:bs-naive-vs-routhian`, Fig. `fig:bs-J-independence`; text "$\Delta z\sim10^{-7}$ ... nearest barrier ... order unity"; identity $F_1(z)/F_3(z)\equiv4$.
- **Status: PASS.** Reproduced $\Delta z \approx 4.49\times10^{-7}$ to $4.52\times10^{-7}$ depending on $r_{\rm UV}$ (Table C); `max|F1/F3-4| = 0.000e+00`.

### Notebook — `notebooks/Notebook_Krylov_Final.ipynb`
- **Role:** exploratory notebook that predates the standalone `quiver_comparison/` scripts (see `run_quiver_comparison.py` docstring: "repeats the Quiver-1 fixed-charge-Routhian-vs-naive analysis (already validated in Notebook_Krylov_Final.ipynb)").
- **Status: PASS (executes without error end-to-end)** via `jupyter nbconvert --execute` on 2026-09-09; not separately cross-checked cell-by-cell against report equations beyond this, since its results are superseded by the audited `quiver_comparison/` scripts above.

### `report/supplementary/Notes.tex`
- **Role:** working derivation notes explicitly cited by code comments (e.g. `blackstring_quiver/run_blackstring_quiver.py`: "the supervisors' instructions (see the shared email quoted in Notes.tex)"; "Notes.tex, Sec. 'Extension: Fixed-charge Routhian for holographic 6d SCFT backgrounds'"). Included for scientific traceability; it is a draft, not the final report.
- **Status: NOT independently verified** (not compiled/checked in this audit — see `REPOSITORY_AUDIT.md`).

---

## Figure provenance summary

| Figure | Type | Reproducing script |
|---|---|---|
| `fig:flowchart` (Krylov construction diagram) | TikZ (in-document) | n/a — drawn directly in LaTeX |
| `fig:toy-ct` (`ct_example.png`) | claimed numerical | **none found** — NOT VERIFIED |
| `fig:ads-cylinder`, `fig:naive-vs-routhian-schematic` | TikZ (in-document) | n/a — drawn directly in LaTeX |
| `fig:ads3-validation` (`holographic_krylov_validation.png`) | numerical | `ads3/ads3_validation.py` — VERIFIED |
| `fig:quiver-f1-comparison`, `fig:quiver-intercept`, `fig:quiver-P-scan` | numerical | `quiver_comparison/run_*.py` — VERIFIED |
| `fig:bs-naive-vs-routhian`, `fig:blackstring-h4`, `fig:bs-J-independence` | numerical | `blackstring_quiver/run_blackstring_quiver.py` — VERIFIED |

No figure in the report is imported from a third-party/cited paper (the report's own diagrams are either TikZ or produced by the scripts above).
