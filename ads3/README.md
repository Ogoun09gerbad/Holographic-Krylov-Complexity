# ads3/ — Global AdS$_3$ analytic benchmark

Corresponds to report section **"Global AdS$_3$: analytic benchmark"** (`sec:ads3`).

## What this computes

`ads3_validation.py` performs two independent checks of the closed-form results derived by hand in the report for a charged point particle in global AdS$_3$ ($f=1+r^2$, $g=1/(1+r^2)$, $h=r^2$):

1. **Symbolic (SymPy).** Re-derives from the Routhian $R=-m_{\mathrm{eff}}(r)\sqrt{f-g\dot r^2}$, $m_{\mathrm{eff}}=\sqrt{m^2+J^2/h}$:
   - the turning-point energy $E^2$ (eq. `ads3-energy-initial`),
   - the first-order radial equation $\dot r^2(r)$ (eq. `ads3-radial-firstorder`),
   - the initial radial acceleration $\ddot r(0)$ (eq. `ads3-acceleration`),
   - the short-time slope $\alpha(J,m,r_0)$ of $P_r^{(R)}(t)$ (eq. `ads3-alpha`),
   - the neutral limit $J\to0$ and the boundary limit $r_0\to\infty$,

   and symbolically compares each to the closed form quoted in the `.tex` source (`sympy.simplify` of the difference).

2. **Numerical (SciPy).** Lambdifies the exact Euler–Lagrange equation from the same Routhian and integrates it with `scipy.integrate.solve_ivp` (method `DOP853`, `rtol=1e-13`, `atol=1e-15`) from a turning point $r(0)=r_0$, $\dot r(0)=0$, for five representative $(m,J,r_0)$. The short-time slope of $P_r^{(R)}(t)$ is extracted by a least-squares fit through the origin on $t\in[10^{-6},2\times10^{-4}]$ and compared to the analytic $\alpha(J,m,r_0)$.

## How to run

```bash
python ads3_validation.py
```

No arguments, no random seed needed (fully deterministic ODE integration). Produces `holographic_krylov_validation.png` in this directory (referenced by the report as `sec:ads3`, `\label{fig:ads3-validation}`) and prints the symbolic PASS/FAIL table and the numeric relative-error table to stdout.

## Verification status

Re-run on 2026-09-09 (Python 3.14.7, numpy 2.5.2, sympy 1.14.0, scipy 1.18.0, matplotlib 3.11.1): all 6 symbolic checks PASS; worst numeric relative error over the 5 test cases was $1.04\times10^{-7}$, consistent with the report's statement "the relative deviation was at most $1.0\times10^{-7}$".
