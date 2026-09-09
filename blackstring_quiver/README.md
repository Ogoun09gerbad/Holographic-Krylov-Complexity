# blackstring_quiver/ — Quadratic black-string-chain background

Corresponds to report section **"Quadratic black-string-chain background"** (`sec:blackstring`).

## Files

| File | Role |
|---|---|
| `geometry.py` | Implements the near-horizon $\mathcal N=(0,4)$ D2–D4–D6–NS5 black-string-chain background of Couzens, Lozano, Petri, Vandoren (arXiv:2109.10413) with the quadratic rank profile $h_4(z)=\alpha+\beta z-\tfrac{\gamma}{2}z^2$, and derives $F_1(z),F_2(z),F_3(z)$ of eq. `blackstring-Fi` (including the exact identity $F_1/F_3\equiv4$). |
| `run_blackstring_quiver.py` | Main driver: well/barrier structure of $h_4(z)$, the naive-vs-Routhian mechanism table, the $J$-independence-of-equilibrium check, and the UV-cutoff-vs-displacement scan. |

`run_blackstring_quiver.py` imports `quiver_comparison/routhian.py` **unchanged** via a relative `sys.path.insert(0, os.path.join(HERE, "..", "quiver_comparison"))`, because the reduced Lagrangian here has exactly the same generic form ($f_1(\eta)\to F_1(z)$ etc.) as the 6d quiver case. **This means `blackstring_quiver/` and `quiver_comparison/` must remain sibling directories** for this script to run — do not move one without the other.

## How to run

```bash
python run_blackstring_quiver.py   # -> figures/01-03*.png, tableA-D_*.csv
```

Requires `../quiver_comparison/routhian.py` to be present (see above).

## Figures actually cited in the report

`01_h4_profile.png`, `02_naive_vs_routhian.png`, `03_J_independent_equilibrium.png` are the three figures referenced by `\IfFileExists` in the report. Two additional PNGs found in the original project directory (`02_J_independent_equilibrium.png`, `03_barrier_crossing_trajectories.png`) were **not** reproduced by the current version of `run_blackstring_quiver.py` when re-run from a clean copy, and are therefore excluded from this repository as stale outputs of an earlier script version — see `REPOSITORY_AUDIT.md`.

## Reference benchmark

$P=6$, $\gamma=32$, $(\beta_0,\ldots,\beta_5)=(22,21,19,17,15,15)$ with $\beta_6$ fixed by the closing condition (Sec. `sec:blackstring`, "Multi-well structure").

## Verification status

Re-run on 2026-09-09: Table A mechanism, Table B $J$-independence, Table C UV-cutoff scan, and Table D minimal-quiver check all reproduced. Displacement scale $\Delta z \approx 4.49\times10^{-7}$ at the reference well (Table B, $J=0$), consistent with the report's "order $10^{-7}$" claim. The identity $F_1/F_3\equiv4$ verified to hold with `max|F1/F3-4| = 0.000e+00`.
