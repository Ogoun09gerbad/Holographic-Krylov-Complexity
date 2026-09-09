"""
The original 'early-time zoom' figure (05_early_time_zoom.png) plots
P_naive(t) and P_R(t) on their natural scale (up to ~40-50 by t=0.4).
On that scale the nonzero naive intercept -- P_naive(0) = |J|/sqrt(f3(eta0)),
only 0.09-0.69 depending on (quiver, J) -- is under 1% of the plot height
and is visually indistinguishable from zero. That is exactly what the
supervisors flagged as "too small to see" on the printed slide.

This script makes the effect actually visible by plotting the DIFFERENCE
    Delta P(t) = P_naive(t) - P_R(t)
which isolates the nonzero-intercept effect from the shared linear growth
common to both momenta. Since P_R(0)=0 exactly, Delta P(0) = P_naive(0) by
construction, so the y-intercept of this plot is precisely the number
quoted in the "Two quivers" table (0.088, 0.263, 0.231, 0.692).

Run with:  python run_intercept_plot.py
Produces:  figures/09_naive_minus_routhian_intercept.png
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import run_quiver_comparison as rqc

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")


def main():
    t_early = np.linspace(1e-4, 0.05, 300)
    out = {}
    for qname, qd in rqc.QUIVERS.items():
        out[qname] = {}
        for J in rqc.J_VALUES:
            sol = rqc.rd.integrate_trajectory(J, rqc.ETA0, rqc.R_UV, 0.05, qd["f1"], qd["f2"], qd["f3"])
            data = rqc.rd.naive_and_routhian_momenta(t_early, sol, J, rqc.ETA0, rqc.R_UV, qd["f1"], qd["f2"], qd["f3"])
            out[qname][J] = data

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=False)
    for ax, (qname, qd) in zip(axes, rqc.QUIVERS.items()):
        colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(rqc.J_VALUES)))
        for J, c in zip(rqc.J_VALUES, colors):
            dP = out[qname][J]["P_naive"] - out[qname][J]["P_R"]
            analytic = abs(J) / np.sqrt(qd["f3"](np.array([rqc.ETA0])))[0]
            ax.plot(t_early, dP, color=c, lw=2, label=f"$J={J:g}$ (analytic: {analytic:.3f})")
            ax.axhline(analytic, color=c, ls=":", lw=1)
        ax.set_title(qd["label"])
        ax.set_xlabel("$t$")
        ax.set_ylabel(r"$P_\rho^{\rm naive}(t)-P_\rho^{(R)}(t)$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    plt.suptitle(r"The naive-momentum intercept, isolated: $\Delta P(t)\to |J|/\sqrt{f_3(\eta_0)}$ as $t\to0$")
    plt.tight_layout()
    outpath = os.path.join(FIGDIR, "09_naive_minus_routhian_intercept.png")
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
