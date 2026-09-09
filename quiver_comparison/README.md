# quiver_comparison/ — Two 6d SCFT quiver geometries

Corresponds to report sections **"Top-down six-dimensional SCFT backgrounds"** (`sec:6d`) and **"Two 6d SCFT quiver geometries"** (`sec:quivers`).

## Files

| File | Role |
|---|---|
| `geometry.py` | Builds the warp factors $f_1(\eta),f_2(\eta),f_3(\eta)$ for Quiver 1 (triangular rank, single flavour node) and Quiver 2 (trapezoidal rank, two flavour nodes) from the rank functions of arXiv:2603.10106, plus a geometry self-validation (continuity of $\alpha(\eta)$ and its derivative at the piecewise breakpoints, positivity checks). |
| `routhian.py` | Generic fixed-charge Routhian reduction for a metric of the form $\mathrm ds^2 = f_1(\eta)\mathrm ds^2_{\mathrm{AdS}_7} + f_2(\eta)\mathrm d\eta^2 + f_3(\eta)\mathrm d\varphi^2$ (Sec. `sec:6d`). Also reused unchanged by `blackstring_quiver/run_blackstring_quiver.py`. |
| `run_quiver_comparison.py` | Main driver: naive-vs-Routhian momentum comparison, parameter scan over $J$, exact $P_r^{(R)}=Ht$ radial check, for both quivers at $N=1,P=100,r_{\rm UV}=\ln(0.01),\eta_0=50$. |
| `run_parameter_family_scan.py` | Scans quiver size $P\in\{5,20,50,100,200\}$ at fixed $N=1$ to test robustness of the fixed-charge mechanism across geometry size (Sec. `sec:quivers`, "Parameter scan"). |
| `run_intercept_plot.py` | Isolates the naive-vs-Routhian intercept effect by plotting $\Delta P(t)=P_{\rm naive}(t)-P_R(t)$, since the raw intercept is visually indistinguishable from zero on the natural plot scale. |
| `sympy_verify_quiver2.py` | Independent symbolic (SymPy) re-derivation of Quiver 2's rank function $\alpha_2(\eta)$ from $\alpha''=-81\pi^2 R_2(\eta)$ plus boundary/continuity conditions, cross-checked against `geometry.py`'s closed form. |

## How to run

```bash
python run_quiver_comparison.py        # -> figures/01-08*.png, results_table.csv, parameter_scan.csv
python run_parameter_family_scan.py    # -> figures/09_f1_shape_vs_P.png, figures/10_scan_vs_P.png, parameter_family_scan.csv
python run_intercept_plot.py           # -> figures/09_naive_minus_routhian_intercept.png (requires figures/ to already exist)
python sympy_verify_quiver2.py         # stdout only, no file output
```

All scripts write output relative to their own file location (`HERE = os.path.dirname(os.path.abspath(__file__))`), so they can be run from any working directory.

## Figures actually cited in the report

Of the 11 figures in `figures/`, three are referenced by `\IfFileExists` in the report: `01_geometry_f1_comparison.png`, `09_naive_minus_routhian_intercept.png`, `10_scan_vs_P.png`. The remaining eight (`02`–`08`, `09_f1_shape_vs_P`) are legitimate reproducible intermediate diagnostics produced by the same scripts but not included in the final report figures.

## Parameter provenance

$N=1$, $P=100$ are taken directly from arXiv:2603.10106. $r_{\rm UV}=\ln(0.01)$ was reconstructed from that reference's figure captions (moderate confidence, as stated in the report). $\eta_0=50$ is the report's own choice, shared between both quivers so the comparison isolates geometry shape rather than starting point.

## Verification status

Re-run on 2026-09-09: naive intercepts $|J|/\sqrt{f_3(\eta_0)}$ reproduced as $0.0877, 0.2630$ (Quiver 1, $J=1,3$) and $0.2307, 0.6920$ (Quiver 2, $J=1,3$) — exact match to the report's Table in `sec:quivers`. The $N$-scaling identity $f_i(\eta;N)=\sqrt N f_i(\eta;1)$ and the exact radial identity $P_r^{(R)}(t)=Ht$ were both confirmed numerically (`max|Pr-Ht|` column in stdout, of order $10^{-12}$).
