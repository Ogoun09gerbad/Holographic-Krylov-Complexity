"""
Top-level driver: repeats the Quiver-1 fixed-charge-Routhian-vs-naive
analysis (already validated in Notebook_Krylov_Final.ipynb) for Quiver 2,
and produces every comparison plot / table requested for the project's
"other quiver" extension.

Run with:  python run_quiver_comparison.py
Produces:  figures/*.png  and prints all numerical tables to stdout
        (also saved to results_table.csv / parameter_scan.csv)
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

import geometry
import routhian as rd

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

# ----------------------------------------------------------------------
# Parameters
#   N, P: TAKEN FROM THE REFERENCE (P=100 appears repeatedly in the
#         figure captions retrieved for both quivers; N is an overall
#         rank scale that does not affect the shape of f1,f2,f3).
#   r_UV: TAKEN FROM THE REFERENCE figure captions (e^{r_UV}=0.01),
#         moderate confidence -- see audit notes in the report.
#   eta0: our own choice (ASSUMED), same value used for both quivers
#         so that the comparison isolates the effect of the geometry
#         rather than of a different starting point.
# ----------------------------------------------------------------------
N_RANK = 1.0
P_QUIVER = 100.0
R_UV = np.log(0.01)          # ~ -4.60517
ETA0 = 50.0
J_VALUES = [0.0, 1.0, 3.0]
T_MAX = 3.0
N_TPLOT = 400

QUIVERS = {
    "quiver1": dict(label="Quiver 1 (triangular, 1 flavour node)"),
    "quiver2": dict(label="Quiver 2 (trapezoidal, 2 flavour nodes)"),
}
for qname, qd in QUIVERS.items():
    f1, f2, f3, f6 = geometry.make_metric_functions(qname, N=N_RANK, P=P_QUIVER)
    qd.update(f1=f1, f2=f2, f3=f3, f6=f6)


def run_all_trajectories():
    t_plot = np.linspace(1e-4, T_MAX, N_TPLOT)
    out = {}
    for qname, qd in QUIVERS.items():
        out[qname] = {}
        for J in J_VALUES:
            sol = rd.integrate_trajectory(J, ETA0, R_UV, T_MAX, qd["f1"], qd["f2"], qd["f3"])
            data = rd.naive_and_routhian_momenta(t_plot, sol, J, ETA0, R_UV, qd["f1"], qd["f2"], qd["f3"])
            data["sol_status"] = sol.status
            out[qname][J] = data
    return t_plot, out


def fit_power_laws(t_plot, P_R, C_R, fit_tmax=0.3):
    mask = t_plot < fit_tmax
    tt = t_plot[mask]

    # P_R(t) = a t + b t^3
    A = np.vstack([tt, tt**3]).T
    coeff_P, *_ = np.linalg.lstsq(A, P_R[mask], rcond=None)

    # C_R(t) = c t^2 + d t^4
    A2 = np.vstack([tt**2, tt**4]).T
    coeff_C, *_ = np.linalg.lstsq(A2, C_R[mask], rcond=None)

    # log-log slope estimate at small t
    logt = np.log(tt[tt > 1e-3])
    logP = np.log(np.abs(P_R[mask][tt > 1e-3]) + 1e-300)
    slope_loglog = np.polyfit(logt, logP, 1)[0]

    return dict(a=coeff_P[0], b=coeff_P[1], c=coeff_C[0], d=coeff_C[1], loglog_slope=slope_loglog)


def main():
    print("=" * 70)
    print("GEOMETRY SELF-VALIDATION")
    print("=" * 70)
    geometry._validate()

    print()
    print("=" * 70)
    print(f"PARAMETERS: N={N_RANK}, P={P_QUIVER}, r_UV={R_UV:.6f}, eta0={ETA0}, "
          f"J in {J_VALUES}, t_max={T_MAX}")
    print("=" * 70)

    t_plot, results = run_all_trajectories()

    # --- Figure 1: geometry comparison f1(eta) ---
    eta_g = np.linspace(1.0, P_QUIVER - 1.0, 600)
    plt.figure(figsize=(8, 5))
    for qname, qd in QUIVERS.items():
        plt.plot(eta_g, qd["f1"](eta_g), label=f"{qd['label']}: $f_1(\\eta)$")
    plt.xlabel(r"$\eta$")
    plt.ylabel(r"$f_1(\eta)$")
    plt.title("Background geometry: $f_1(\\eta)$ for both quivers")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "01_geometry_f1_comparison.png"), dpi=150)
    plt.close()

    plt.figure(figsize=(12, 4))
    for i, (qname, qd) in enumerate(QUIVERS.items()):
        plt.subplot(1, 2, i + 1)
        plt.plot(eta_g, qd["f2"](eta_g), label="$f_2$")
        plt.plot(eta_g, qd["f3"](eta_g), label="$f_3$")
        plt.xlabel(r"$\eta$")
        plt.title(qd["label"])
        plt.legend()
        plt.grid(alpha=0.3)
    plt.suptitle("$f_2(\\eta), f_3(\\eta)$ for both quivers")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "02_geometry_f2_f3_comparison.png"), dpi=150)
    plt.close()

    # --- Figure 2: eta(t) trajectories ---
    plt.figure(figsize=(12, 5))
    for i, (qname, qd) in enumerate(QUIVERS.items()):
        plt.subplot(1, 2, i + 1)
        for J in J_VALUES:
            plt.plot(t_plot, results[qname][J]["eta"], label=f"$J={J}$")
        plt.axhline(ETA0, color="k", ls=":", lw=1)
        plt.xlabel("$t$")
        plt.ylabel(r"$\eta(t)$")
        plt.title(qd["label"])
        plt.legend()
        plt.grid(alpha=0.3)
    plt.suptitle(r"Quiver-direction trajectory $\eta(t)$ for several $R$-charges")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "03_eta_trajectories.png"), dpi=150)
    plt.close()

    # --- Figure 3: naive vs Routhian, full range ---
    plt.figure(figsize=(12, 5))
    for i, (qname, qd) in enumerate(QUIVERS.items()):
        plt.subplot(1, 2, i + 1)
        colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(J_VALUES)))
        for J, c in zip(J_VALUES, colors):
            plt.plot(t_plot, results[qname][J]["P_naive"], color=c, label=f"naive, $J={J}$")
            plt.plot(t_plot, results[qname][J]["P_R"], color=c, ls="--", label=f"Routhian, $J={J}$")
        plt.xlabel("$t$")
        plt.ylabel(r"$P_\rho(t)$")
        plt.title(qd["label"])
        plt.legend(fontsize=7, ncol=2)
        plt.grid(alpha=0.3)
    plt.suptitle("Naive vs. fixed-charge (Routhian) proper momentum")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "04_naive_vs_routhian.png"), dpi=150)
    plt.close()

    # --- Figure 4: early-time zoom ---
    mask_early = t_plot < 0.4
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for i, (qname, qd) in enumerate(QUIVERS.items()):
        colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(J_VALUES)))
        for J, c in zip(J_VALUES, colors):
            axes[i, 0].plot(t_plot[mask_early], results[qname][J]["P_naive"][mask_early], color=c, label=f"$J={J}$")
            axes[i, 1].plot(t_plot[mask_early], results[qname][J]["P_R"][mask_early], color=c, label=f"$J={J}$")
        axes[i, 0].set_title(f"{qd['label']}: naive $P_\\rho^{{\\rm naive}}(t)$")
        axes[i, 1].set_title(f"{qd['label']}: fixed-charge $P_\\rho^{{(R)}}(t)$")
        for ax in axes[i]:
            ax.set_xlabel("$t$")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=8)
    plt.suptitle("Early-time behaviour: nonzero naive intercept vs. Routhian linear vanishing")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "05_early_time_zoom.png"), dpi=150)
    plt.close()

    # --- Figure 5: log-log scaling + Figure 6: C(t) ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fit_results = {}
    for i, (qname, qd) in enumerate(QUIVERS.items()):
        colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(J_VALUES)))
        fit_results[qname] = {}
        for J, c in zip(J_VALUES, colors):
            P_R = np.abs(results[qname][J]["P_R"])
            C_R = cumulative_trapezoid(P_R, t_plot, initial=0)
            fits = fit_power_laws(t_plot, results[qname][J]["P_R"], C_R)
            fit_results[qname][J] = fits

            axes[i, 0].loglog(t_plot[t_plot > 1e-3], P_R[t_plot > 1e-3], color=c, label=f"$J={J}$ (slope={fits['loglog_slope']:.3f})")
            axes[i, 1].plot(t_plot, C_R, color=c, label=f"$J={J}$")
        axes[i, 0].set_title(f"{qd['label']}: $|P_R(t)|$ log-log")
        axes[i, 0].set_xlabel("$t$")
        axes[i, 0].set_ylabel(r"$|P_R(t)|$")
        axes[i, 0].legend(fontsize=8)
        axes[i, 0].grid(alpha=0.3, which="both")
        axes[i, 1].set_title(f"{qd['label']}: $C_R(t)=\\int_0^t |P_R|dt'$")
        axes[i, 1].set_xlabel("$t$")
        axes[i, 1].set_ylabel(r"$C_R(t)$")
        axes[i, 1].legend(fontsize=8)
        axes[i, 1].grid(alpha=0.3)
    plt.suptitle(r"Power-law check: $P_R(t)\sim t^1$, $C_R(t)\sim t^2$")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "06_loglog_and_complexity.png"), dpi=150)
    plt.close()

    # --- Figure 7: exact P_r(t) = H t check ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for i, (qname, qd) in enumerate(QUIVERS.items()):
        colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(J_VALUES)))
        for J, c in zip(J_VALUES, colors):
            Pr_t = results[qname][J]["P_r"]
            H = results[qname][J]["H"]
            axes[i].plot(t_plot, Pr_t - H * t_plot, color=c, label=f"$J={J}$")
        axes[i].set_title(f"{qd['label']}: residual $P_r(t)-Ht$")
        axes[i].set_xlabel("$t$")
        axes[i].grid(alpha=0.3)
        axes[i].legend(fontsize=8)
    plt.suptitle("Exact identity check: $P_r(t)=Ht$ (should be ~0 to numerical precision)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "07_exact_Pr_Ht_check.png"), dpi=150)
    plt.close()

    # --- Numerical table: intercepts, slopes, comparison with analytics ---
    print()
    print("=" * 100)
    print("TABLE: naive/Routhian intercepts and slopes, numeric vs analytic (Notes.tex eqs. 6.21-6.24)")
    print("=" * 100)
    rows = []
    header = (f"{'quiver':>8} {'J':>5} {'P_naive(0) num':>16} {'|J|/sqrt(f3)':>14} "
              f"{'P_R(0) num':>12} {'P_R slope num':>15} {'A_6d analytic':>15} "
              f"{'max|Pr-Ht|':>12}")
    print(header)
    for qname, qd in QUIVERS.items():
        for J in J_VALUES:
            d = results[qname][J]
            short = rd.short_time_coefficients(J, ETA0, R_UV, qd["f1"], qd["f2"], qd["f3"])
            slope_num = (d["P_R"][10] - d["P_R"][0]) / (t_plot[10] - t_plot[0])
            max_resid = np.max(np.abs(d["P_r"] - d["H"] * t_plot))
            row = dict(quiver=qname, J=J,
                       P_naive_0=d["P_naive"][0], naive_intercept_analytic=short["naive_intercept"],
                       P_R_0=d["P_R"][0], slope_num=slope_num, A_6d=abs(short["A_6d"]),
                       max_resid=max_resid, eta_min=d["eta"].min(), eta_max=d["eta"].max())
            rows.append(row)
            print(f"{qname:>8} {J:5.1f} {d['P_naive'][0]:16.6f} {short['naive_intercept']:14.6f} "
                  f"{d['P_R'][0]:12.6f} {slope_num:15.6f} {abs(short['A_6d']):15.6f} {max_resid:12.3e}")

    with open(os.path.join(HERE, "results_table.csv"), "w", newline="") as fcsv:
        writer = csv.DictWriter(fcsv, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # --- Parameter scan over a finer grid of J, for both quivers ---
    print()
    print("=" * 100)
    print("PARAMETER SCAN over J (falling / eta-confinement, intercepts, slopes, power law)")
    print("=" * 100)
    J_scan = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0]
    scan_rows = []
    scan_header = (f"{'quiver':>8} {'J':>6} {'eta_min':>9} {'eta_max':>9} {'confined?':>10} "
                   f"{'P_naive(0)':>12} {'P_R(0)':>10} {'P_R slope':>11} {'loglog slope':>13} "
                   f"{'b_intercept(fit)':>18}")
    print(scan_header)
    for qname, qd in QUIVERS.items():
        for J in J_scan:
            sol = rd.integrate_trajectory(J, ETA0, R_UV, T_MAX, qd["f1"], qd["f2"], qd["f3"])
            d = rd.naive_and_routhian_momenta(t_plot, sol, J, ETA0, R_UV, qd["f1"], qd["f2"], qd["f3"])
            short = rd.short_time_coefficients(J, ETA0, R_UV, qd["f1"], qd["f2"], qd["f3"])
            slope_num = (d["P_R"][10] - d["P_R"][0]) / (t_plot[10] - t_plot[0])

            P_R = np.abs(d["P_R"])
            C_R = cumulative_trapezoid(P_R, t_plot, initial=0)
            fits = fit_power_laws(t_plot, d["P_R"], C_R)

            confined = (d["eta"].min() > 0.5) and (d["eta"].max() < P_QUIVER - 0.5)
            row = dict(quiver=qname, J=J, eta_min=d["eta"].min(), eta_max=d["eta"].max(),
                       confined=confined, P_naive_0=d["P_naive"][0], P_R_0=d["P_R"][0],
                       slope_num=slope_num, loglog_slope=fits["loglog_slope"], b_fit=fits["b"])
            scan_rows.append(row)
            print(f"{qname:>8} {J:6.1f} {d['eta'].min():9.3f} {d['eta'].max():9.3f} {str(confined):>10} "
                  f"{d['P_naive'][0]:12.6f} {d['P_R'][0]:10.6f} {slope_num:11.6f} "
                  f"{fits['loglog_slope']:13.4f} {fits['b']:18.4e}")

    with open(os.path.join(HERE, "parameter_scan.csv"), "w", newline="") as fcsv:
        writer = csv.DictWriter(fcsv, fieldnames=list(scan_rows[0].keys()))
        writer.writeheader()
        writer.writerows(scan_rows)

    # --- Figure 8: parameter scan plot ---
    plt.figure(figsize=(9, 6))
    for qname, qd in QUIVERS.items():
        Js = [r["J"] for r in scan_rows if r["quiver"] == qname]
        p0s = [r["P_naive_0"] for r in scan_rows if r["quiver"] == qname]
        pR0s = [r["P_R_0"] for r in scan_rows if r["quiver"] == qname]
        plt.plot(Js, p0s, "o-", label=f"{QUIVERS[qname]['label']}: $P_{{naive}}(0)$")
        plt.plot(Js, pR0s, "s--", label=f"{QUIVERS[qname]['label']}: $P_R(0)$")
    plt.xlabel("$J$")
    plt.ylabel(r"intercept at $t=0$")
    plt.title(r"Parameter scan: naive intercept grows with $J$, Routhian intercept stays $\approx 0$")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "08_parameter_scan.png"), dpi=150)
    plt.close()

    print()
    print(f"All figures written to: {FIGDIR}")
    print(f"Tables written to: results_table.csv, parameter_scan.csv (in {HERE})")


if __name__ == "__main__":
    main()
