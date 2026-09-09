"""
Parameter-family scan: within Quiver 1 and Quiver 2 (arXiv:2603.10106,
Sec. 2.1.1-2.1.2 -- the ONLY two quivers this reference actually defines),
vary the number of gauge nodes P to obtain several genuinely different
background geometries per quiver family. N (overall rank scale) is NOT
scanned: geometry.py's own self-check (see run_quiver_comparison.py output
and the N-scaling check below) shows f1,f2,f3 all scale as exactly sqrt(N)
for both quivers -- a pure normalization, not a shape change -- so varying
N alone would not produce a new geometry, only rescale the whole problem.

This is NOT an attempt to reproduce a "Quiver 3" (the reference does not
have one); it is a systematic robustness check of the fixed-charge
Routhian mechanism across quiver SIZE, complementing the quiver TYPE
comparison of run_quiver_comparison.py.
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import geometry
import routhian as rd

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

N_RANK = 1.0
P_VALUES = [5.0, 20.0, 50.0, 100.0, 200.0]
R_UV = np.log(0.01)
J_VALUES = [0.0, 1.0, 3.0]
T_MAX = 3.0
N_TPLOT = 400


def eta0_of_P(P):
    """Our own choice (ASSUMED): probe released at the midpoint of the
    quiver domain [0, P+1], so the comparison automatically adapts to
    each quiver's own size instead of using a fixed eta0 that would sit
    outside the domain for small P."""
    return (P + 1.0) / 2.0


def verify_N_scaling():
    print("=== N-scaling check (analytic claim: f1,f2,f3 ~ sqrt(N), shape-invariant) ===")
    eta = np.array([30.0])
    for qname in ["quiver1", "quiver2"]:
        f1a, f2a, f3a, _ = geometry.make_metric_functions(qname, N=1.0, P=100.0)
        f1b, f2b, f3b, _ = geometry.make_metric_functions(qname, N=4.0, P=100.0)
        print(f"  {qname}: f1(N=4)/f1(N=1)={f1b(eta)[0]/f1a(eta)[0]:.6f}, "
              f"f2 ratio={f2b(eta)[0]/f2a(eta)[0]:.6f}, f3 ratio={f3b(eta)[0]/f3a(eta)[0]:.6f} "
              f"(all should equal sqrt(4)=2.000000)")
    print()


def main():
    verify_N_scaling()

    # --- Figure: f1(eta) shape across P, both quivers, normalized xi=eta/(P+1) ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for i, qname in enumerate(["quiver1", "quiver2"]):
        colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(P_VALUES)))
        for P, c in zip(P_VALUES, colors):
            f1, f2, f3, f6 = geometry.make_metric_functions(qname, N=N_RANK, P=P)
            eta_g = np.linspace(0.02 * (P + 1), 0.98 * (P + 1), 400)
            xi = eta_g / (P + 1)
            axes[i].plot(xi, f1(eta_g) / np.sqrt(N_RANK), color=c, label=f"$P={P:.0f}$")
        axes[i].set_xlabel(r"$\xi=\eta/(P+1)$")
        axes[i].set_ylabel(r"$f_1(\eta)/\sqrt{N}$")
        axes[i].set_title(f"{qname}: shape of $f_1$ across quiver size $P$")
        axes[i].legend(fontsize=8)
        axes[i].grid(alpha=0.3)
    plt.suptitle("Geometry robustness: $f_1(\\eta)$ shape for 5 quiver sizes per family ($N=1$)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "09_f1_shape_vs_P.png"), dpi=150)
    plt.close()

    # --- Full naive-vs-Routhian check for every (quiver, P, J) combination ---
    print("=" * 110)
    print("PARAMETER-FAMILY SCAN: naive vs Routhian intercepts/slopes across quiver SIZE (P)")
    print("=" * 110)
    header = (f"{'quiver':>8} {'P':>6} {'eta0':>7} {'J':>5} {'eta_max/domain':>15} "
              f"{'P_naive(0)':>12} {'analytic':>10} {'P_R(0)':>10} {'slope num':>11} "
              f"{'|A_6d|':>10} {'max|Pr-Ht|':>12}")
    print(header)

    rows = []
    t_plot = np.linspace(1e-4, T_MAX, N_TPLOT)

    for qname in ["quiver1", "quiver2"]:
        for P in P_VALUES:
            f1, f2, f3, f6 = geometry.make_metric_functions(qname, N=N_RANK, P=P)
            eta0 = eta0_of_P(P)

            eta_scan = np.linspace(1.0, P - 1.0 if P > 2 else P + 0.5, 300)
            eta_scan = eta_scan[(eta_scan > 0) & (eta_scan < P + 1)]
            f1_scan = f1(eta_scan)
            eta_max = eta_scan[np.argmax(f1_scan)]

            for J in J_VALUES:
                sol = rd.integrate_trajectory(J, eta0, R_UV, T_MAX, f1, f2, f3, max_step=0.01)
                d = rd.naive_and_routhian_momenta(t_plot, sol, J, eta0, R_UV, f1, f2, f3)
                short = rd.short_time_coefficients(J, eta0, R_UV, f1, f2, f3)
                slope_num = (d["P_R"][10] - d["P_R"][0]) / (t_plot[10] - t_plot[0])
                max_resid = np.max(np.abs(d["P_r"] - d["H"] * t_plot))

                row = dict(quiver=qname, P=P, eta0=eta0, J=J,
                           eta_max_frac=eta_max / (P + 1),
                           P_naive_0=d["P_naive"][0], naive_analytic=short["naive_intercept"],
                           P_R_0=d["P_R"][0], slope_num=slope_num, A_6d=abs(short["A_6d"]),
                           max_resid=max_resid)
                rows.append(row)
                print(f"{qname:>8} {P:6.0f} {eta0:7.1f} {J:5.1f} {eta_max/(P+1):15.4f} "
                      f"{d['P_naive'][0]:12.6f} {short['naive_intercept']:10.6f} "
                      f"{d['P_R'][0]:10.6f} {slope_num:11.4f} {abs(short['A_6d']):10.4f} "
                      f"{max_resid:12.3e}")

    with open(os.path.join(HERE, "parameter_family_scan.csv"), "w", newline="") as fcsv:
        writer = csv.DictWriter(fcsv, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # --- Summary figure: naive intercept and Routhian slope vs P, for J=3 ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for qname, marker in zip(["quiver1", "quiver2"], ["o-", "s--"]):
        Ps = [r["P"] for r in rows if r["quiver"] == qname and r["J"] == 3.0]
        naive0 = [r["P_naive_0"] for r in rows if r["quiver"] == qname and r["J"] == 3.0]
        PR0 = [r["P_R_0"] for r in rows if r["quiver"] == qname and r["J"] == 3.0]
        slopes = [r["slope_num"] for r in rows if r["quiver"] == qname and r["J"] == 3.0]
        A6ds = [r["A_6d"] for r in rows if r["quiver"] == qname and r["J"] == 3.0]
        alt_marker = marker[0] + ":"
        axes[0].plot(Ps, naive0, marker, label=f"{qname}: $P^{{\\rm naive}}(0)$")
        axes[0].plot(Ps, PR0, alt_marker, label=f"{qname}: $P^{{(R)}}(0)$")
        axes[1].plot(Ps, slopes, marker, label=f"{qname}: numeric slope")
        axes[1].plot(Ps, A6ds, alt_marker, label=f"{qname}: $|A_{{6d}}|$ analytic")
    axes[0].set_xlabel("$P$ (quiver size)")
    axes[0].set_ylabel("intercept at $t=0$")
    axes[0].set_title(r"Naive vs. Routhian intercept, $J=3$, across quiver size")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.3)
    axes[1].set_xlabel("$P$ (quiver size)")
    axes[1].set_ylabel(r"slope of $P_\rho^{(R)}(t)$")
    axes[1].set_title(r"Routhian slope: numeric vs. analytic $A_{6d}$, across quiver size")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.3)
    plt.suptitle("Robustness of the fixed-charge mechanism across quiver size $P$ ($J=3$)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "10_scan_vs_P.png"), dpi=150)
    plt.close()

    print()
    print(f"Figures written to: {FIGDIR} (09_f1_shape_vs_P.png, 10_scan_vs_P.png)")
    print(f"Table written to: {os.path.join(HERE, 'parameter_family_scan.csv')}")


if __name__ == "__main__":
    main()
