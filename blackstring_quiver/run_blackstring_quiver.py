"""
Fixed-charge Routhian analysis for the QUADRATIC black-string-chain quiver
of arXiv:2109.10413 (Couzens, Lozano, Petri, Vandoren), following exactly
the supervisors' instructions (see the shared email quoted in Notes.tex):

    "I would make the point particle move along the S^2, the t and r of the
    AdS3 and the z coordinate (eta in the arXiv:2512.14812 paper). You
    should use the Routhian as you have and see whether there is
    interesting behaviour because of this quadratic term in the rank
    function."

Because the fixed-charge Routhian machinery of quiver_comparison/routhian.py
was derived (Notes.tex, Sec. "Extension: Fixed-charge Routhian for
holographic 6d SCFT backgrounds") for a GENERIC metric of the form
    ds_E^2 = f1(eta)ds^2_AdS + f2(eta)deta^2 + f3(eta)dphi^2 (+ frozen theta),
and this new background has exactly that structure (see geometry.py's
docstring for the derivation of F1,F2,F3 from the arXiv:2109.10413 metric),
we reuse quiver_comparison/routhian.py UNCHANGED, simply feeding it the new
F1,F2,F3 of this quiver instead of Quiver 1/Quiver 2's f1,f2,f3.

DESIGN NOTE (learned from an earlier, much slower version of this script):
a naive grid scan over (z0, r_UV, J) with a long, uniformly-sampled t-grid
is wasteful, because (a) each trajectory only needs to be integrated ONCE
(figures should reuse the same `sol` used for the tables, not re-integrate),
and (b) the "kick" magnitude H(J,eta0,r_UV) that drives the eta-motion sets
its OWN internal relaxation timescale ~1/H; for the reference quiver
H is O(1-1e4) already at r_UV=ln(0.01), so a t-grid resolved at t~0.1-1 is
already deep in the "frozen" late-time regime -- there is no need for
t_max > a few units, and no need for J or r_UV scans finer than a handful
of representative points to see the qualitative (saturating) behaviour.

Run with: python run_blackstring_quiver.py
"""

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema

import geometry as bg

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

sys.path.insert(0, os.path.join(HERE, "..", "quiver_comparison"))
import routhian as rd  # noqa: E402  (generic in f1,f2,f3; TAKEN FROM quiver_comparison/, unmodified)

R_UV = np.log(0.01)
J_VALUES = [0.0, 1.0, 3.0]
T_MAX = 6.0
N_TPLOT = 600


def find_wells_and_barriers(q, n=40001):
    z = np.linspace(1e-4, q.z_max - 1e-4, n)
    h = q.h4(z)
    maxima = argrelextrema(h, np.greater)[0]
    minima = argrelextrema(h, np.less)[0]
    wells = [(z[i], h[i], int(q.segment_index(np.array([z[i]]))[0])) for i in maxima]
    barriers = [(z[i], h[i], int(q.segment_index(np.array([z[i]]))[0])) for i in minima]
    return wells, barriers, z, h


def main():
    print("=" * 70)
    print("GEOMETRY SELF-VALIDATION")
    print("=" * 70)
    bg._validate()
    q = bg.make_reference_quiver()
    F1, F2, F3, edPhi = q.metric_functions()

    wells, barriers, z_g, h4_g = find_wells_and_barriers(q)
    print()
    print("Wells (local maxima of h4 == of the J-independent eta-equilibria):")
    for zw, hw, kw in wells:
        print(f"   segment k={kw}: z*={zw:.4f}, h4(z*)={hw:.3f}")
    print("Barriers (local minima of h4, at the interior flavour-brane kinks):")
    for zb, hb, kb in barriers:
        print(f"   kink k={kb}: z_kink={zb:.4f}, h4(z_kink)={hb:.3f}")

    # --- Figure 1: h4(z) with wells/barriers marked ---
    plt.figure(figsize=(10, 5))
    plt.plot(z_g, h4_g, color="C0", lw=1.5, label=r"$h_4(z)$")
    for k in range(1, q.P + 1):
        plt.axvline(2 * np.pi * k, color="gray", ls=":", lw=0.8)
    for zw, hw, kw in wells:
        plt.plot(zw, hw, "o", color="C2")
    for zb, hb, kb in barriers:
        plt.plot(zb, hb, "x", color="C3")
    plt.xlabel(r"$z$")
    plt.ylabel(r"$h_4(z)$")
    plt.title(r"Quadratic rank function $h_4(z)$ (arXiv:2109.10413 Fig.1 quiver): "
              r"wells (o) and kink-barriers (x)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "01_h4_profile.png"), dpi=150)
    plt.close()

    # --- Part A: naive-vs-Routhian mechanism table at the segment-3 well ---
    z0_ref = wells[3][0]
    print()
    print("=" * 100)
    print(f"TABLE A: naive vs Routhian mechanism at z0={z0_ref:.4f} (segment k=3 well), r_UV={R_UV:.6f}")
    print("=" * 100)
    t_plot = np.linspace(1e-4, T_MAX, N_TPLOT)
    header = (f"{'J':>5} {'P_naive(0)':>12} {'|J|/sqrt(F3)':>14} {'P_R(0)':>10} "
              f"{'slope num':>11} {'|A_6d|':>10} {'max|Pr-Ht|':>12}")
    print(header)
    rowsA = []
    trajA = {}
    for J in J_VALUES:
        sol = rd.integrate_trajectory(J, z0_ref, R_UV, T_MAX, F1, F2, F3, max_step=0.01)
        trajA[J] = sol
        d = rd.naive_and_routhian_momenta(t_plot, sol, J, z0_ref, R_UV, F1, F2, F3)
        short = rd.short_time_coefficients(J, z0_ref, R_UV, F1, F2, F3)
        slope_num = (d["P_R"][10] - d["P_R"][0]) / (t_plot[10] - t_plot[0])
        max_resid = np.max(np.abs(d["P_r"] - d["H"] * t_plot))
        rowsA.append(dict(J=J, P_naive_0=d["P_naive"][0], naive_analytic=short["naive_intercept"],
                           P_R_0=d["P_R"][0], slope_num=slope_num, A_6d=abs(short["A_6d"]),
                           max_resid=max_resid))
        print(f"{J:5.1f} {d['P_naive'][0]:12.6f} {short['naive_intercept']:14.6f} "
              f"{d['P_R'][0]:10.6f} {slope_num:11.6f} {abs(short['A_6d']):10.6f} {max_resid:12.3e}")
    with open(os.path.join(HERE, "tableA_mechanism.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rowsA[0].keys()))
        w.writeheader()
        w.writerows(rowsA)

    # --- Figure 2: naive vs Routhian momentum, early-time zoom (reuses trajA) ---
    plt.figure(figsize=(7, 5))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(J_VALUES)))
    mask_early = t_plot < 0.05
    for J, c in zip(J_VALUES, colors):
        d = rd.naive_and_routhian_momenta(t_plot, trajA[J], J, z0_ref, R_UV, F1, F2, F3)
        plt.plot(t_plot[mask_early], d["P_naive"][mask_early], color=c, label=f"naive, $J={J}$")
        plt.plot(t_plot[mask_early], d["P_R"][mask_early], color=c, ls="--", label=f"Routhian, $J={J}$")
    plt.xlabel("$t$")
    plt.ylabel(r"$P_\rho(t)$")
    plt.title("Black-string-chain quiver: naive vs. fixed-charge proper momentum (early time)")
    plt.legend(fontsize=7, ncol=2)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "02_naive_vs_routhian.png"), dpi=150)
    plt.close()

    # --- Part B: J-independence of the eta-equilibrium (F1/F3=4 identically) ---
    print()
    print("=" * 100)
    print("TABLE B: J-independence of the within-well equilibrium (released 1 unit off-centre)")
    print("=" * 100)
    z0_offcentre = z0_ref - 1.0
    T_LONG = 3.0
    t_long = np.linspace(1e-6, T_LONG, 2000)
    print(f"z0 = {z0_offcentre:.4f} (well maximum z*={z0_ref:.4f}), T_max={T_LONG}")
    print(f"{'J':>5} {'z(T_max)':>14} {'total disp. Delta z':>20} {'still in segment 3?':>20}")
    rowsB = []
    trajB = {}
    for J in [0.0, 1.0, 3.0, 8.0, 15.0]:
        sol = rd.integrate_trajectory(J, z0_offcentre, R_UV, T_LONG, F1, F2, F3, max_step=0.002)
        trajB[J] = sol
        z_end = float(sol.sol(np.array([T_LONG]))[0][0])
        disp = z_end - z0_offcentre
        seg = int(q.segment_index(np.array([z_end]))[0])
        in3 = (2 * np.pi * 3 < z_end < 2 * np.pi * 4)
        rowsB.append(dict(J=J, z_end=z_end, delta_z=disp, segment=seg, confined=in3))
        print(f"{J:5.1f} {z_end:14.10f} {disp:20.6e} {str(in3):>20}")
    with open(os.path.join(HERE, "tableB_J_independence.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rowsB[0].keys()))
        w.writeheader()
        w.writerows(rowsB)

    # --- Figure 3: z(t)-z0 (displacement) for table B, log scale to show saturation ---
    plt.figure(figsize=(8, 5))
    for J in [0.0, 1.0, 3.0, 8.0, 15.0]:
        zt = trajB[J].sol(t_long)[0]
        plt.plot(t_long, zt - z0_offcentre, label=f"$J={J}$")
    plt.xlabel("$t$")
    plt.ylabel(r"$z(t)-z_0$")
    plt.title(r"Displacement toward the SAME equilibrium $z^*$ for every $J$"
              "\n" r"(exact consequence of $F_1(z)/F_3(z)\equiv4$ for this background)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "03_J_independent_equilibrium.png"), dpi=150)
    plt.close()

    # --- Part C: does H (the "kick" strength) change the TOTAL displacement,
    #     or only the TIME over which it is reached? (diagnostic run at J=0) ---
    print()
    print("=" * 100)
    print("TABLE C: total displacement Delta z as a function of the UV cutoff r_UV (J=0)")
    print("=" * 100)
    print(f"{'r_UV':>10} {'H':>12} {'Delta z (t=2)':>16}")
    rowsC = []
    for rUVc in [np.log(0.01), -9.2, -20.0]:
        sol = rd.integrate_trajectory(0.0, z0_offcentre, rUVc, 2.0, F1, F2, F3, max_step=0.001)
        H = rd.H_of_J(0.0, z0_offcentre, rUVc, F1, F3)
        z_end = float(sol.sol(np.array([2.0]))[0][0])
        disp = z_end - z0_offcentre
        rowsC.append(dict(r_UV=rUVc, H=H, delta_z=disp))
        print(f"{rUVc:10.3f} {H:12.4e} {disp:16.6e}")
    with open(os.path.join(HERE, "tableC_kick_vs_displacement.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rowsC[0].keys()))
        w.writeheader()
        w.writerows(rowsC)

    # --- Part D: does a SMALLER (minimal-charge) quiver relax the damping? ---
    print()
    print("=" * 100)
    print("TABLE D: minimal-charge quiver (P=1, gamma=2) -- same diagnostic")
    print("=" * 100)
    q_small = bg.BlackStringQuiver(P=1, gamma=2, beta_free=[2], label="minimal P=1,gamma=2")
    F1s, F2s, F3s, _ = q_small.metric_functions()
    wells_s, barriers_s, zg_s, h4g_s = find_wells_and_barriers(q_small, n=4001)
    print("wells:", wells_s, "barriers:", barriers_s)
    z0_ref_s = wells_s[0][0]
    z0_off_s = z0_ref_s - 0.5
    rowsD = []
    for rUVc in [np.log(0.01), -9.2]:
        sol = rd.integrate_trajectory(0.0, z0_off_s, rUVc, 2.0, F1s, F2s, F3s, max_step=0.001)
        H = rd.H_of_J(0.0, z0_off_s, rUVc, F1s, F3s)
        z_end = float(sol.sol(np.array([2.0]))[0][0])
        disp = z_end - z0_off_s
        rowsD.append(dict(quiver="minimal_P1_gamma2", r_UV=rUVc, H=H, delta_z=disp,
                           F2_at_well=float(F2s(np.array([z0_ref_s]))[0])))
        print(f"r_UV={rUVc:10.3f} H={H:12.4e} Delta_z={disp:16.6e} "
              f"F2(z*)={float(F2s(np.array([z0_ref_s]))[0]):.3f}")
    with open(os.path.join(HERE, "tableD_minimal_quiver.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rowsD[0].keys()))
        w.writeheader()
        w.writerows(rowsD)

    print()
    print(f"Figures written to: {FIGDIR}")
    print(f"Tables written to: tableA_mechanism.csv, tableB_J_independence.csv, "
          f"tableC_kick_vs_displacement.csv, tableD_minimal_quiver.csv (in {HERE})")


if __name__ == "__main__":
    main()
