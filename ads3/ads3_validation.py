"""
Independent symbolic + numerical validation of the global-AdS3 benchmark
of the report (Sec. "Global AdS3: analytic benchmark").

Two independent checks are performed:

1. SYMBOLIC (SymPy): the boxed formulas of the report --
   E^2 at the turning point (eq. ads3-energy-initial), the first-order
   radial equation (eq. ads3-radial-firstorder), the initial radial
   acceleration ddot r(0) (eq. ads3-acceleration), and the short-time
   coefficient alpha = A(J,m,r0) of P_r^(R)(t) (eq. ads3-alpha) -- are
   re-derived from scratch from the Routhian
       R(r,rdot;J) = -meff(r) sqrt(f(r) - g(r) rdot^2),
       meff(r) = sqrt(m^2 + J^2/h(r)),
   with f=1+r^2, g=1/(1+r^2), h=r^2, and compared symbolically
   (sympy.simplify of the difference) against the closed forms quoted
   in the .tex source. This is a proof-level check, not a numerical fit.

2. NUMERICAL (SciPy): the exact Euler-Lagrange equation of motion
   rddot(r,rdot;m,J) is derived symbolically (same Routhian) and
   lambdified, then integrated with an adaptive high-order Runge-Kutta
   method (DOP853, rtol=atol~1e-13) from the turning point r(0)=r0,
   rdot(0)=0. The short-time slope of P_r^(R)(t) is extracted by a
   least-squares fit through the origin on an early time window and
   compared with the analytic alpha(J,m,r0) of eq. (ads3-alpha).

Run:  python ads3_validation.py
Produces: holographic_krylov_validation.png (regenerated in place) and
prints a table of relative errors between the numerical slope and the
analytic short-time coefficient, which is the number quoted in the report.
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# 1. Symbolic re-derivation and cross-check against the report formulas
# ----------------------------------------------------------------------

def symbolic_checks():
    print("=== Symbolic cross-checks (SymPy) ===")
    r, rdot, m, J, r0 = sp.symbols('r rdot m J r0', positive=True, real=True)

    f = 1 + r**2
    g = 1 / (1 + r**2)
    h = r**2
    meff = sp.sqrt(m**2 + J**2 / h)
    R = -meff * sp.sqrt(f - g * rdot**2)

    Pr = sp.diff(R, rdot)
    Energy = sp.simplify(rdot * Pr - R)

    # eq. (ads3-energy-initial): E^2 at the turning point r=r0, rdot=0
    E2_0_derived = sp.simplify(Energy.subs({rdot: 0, r: r0})**2)
    E2_0_claimed = (m**2 + J**2 / r0**2) * (1 + r0**2)
    ok1 = sp.simplify(E2_0_derived - E2_0_claimed) == 0
    print(f"[{'PASS' if ok1 else 'FAIL'}] E^2 at turning point matches eq. (ads3-energy-initial)")

    # eq. (ads3-radial-firstorder): rdot^2(r) at fixed E
    E = sp.symbols('E', positive=True)
    sol = sp.solve(sp.Eq(Energy**2, E**2), rdot**2)
    rdot2_derived = sp.simplify(sol[0])
    rdot2_claimed = (1 + r**2)**2 * (1 - (m**2 + J**2 / r**2) * (1 + r**2) / E**2)
    ok2 = sp.simplify(rdot2_derived - rdot2_claimed) == 0
    print(f"[{'PASS' if ok2 else 'FAIL'}] rdot^2(r) matches eq. (ads3-radial-firstorder)")

    # eq. (ads3-acceleration): ddot r(0) = (1/2) d(rdot^2)/dr at r=r0, E=E2_0
    F = rdot2_claimed.subs(E, sp.sqrt(E2_0_claimed))
    rddot0_derived = sp.simplify(sp.diff(F, r).subs(r, r0) / 2)
    rddot0_claimed = (1 + r0**2) * (J**2 - m**2 * r0**4) / (r0 * (m**2 * r0**2 + J**2))
    ok3 = sp.simplify(rddot0_derived - rddot0_claimed) == 0
    print(f"[{'PASS' if ok3 else 'FAIL'}] ddot r(0) matches eq. (ads3-acceleration)")

    # eq. (ads3-alpha): short-time slope alpha of P_r^(R)(t) = alpha t + O(t^3)
    meff0 = sp.sqrt(m**2 + J**2 / r0**2)
    g0 = 1 / (1 + r0**2)
    f0 = 1 + r0**2
    alpha_derived = sp.simplify(meff0 * g0 * rddot0_claimed / sp.sqrt(f0))
    alpha_claimed = (J**2 - m**2 * r0**4) / (r0**2 * sp.sqrt((m**2 * r0**2 + J**2) * (1 + r0**2)))
    ok4 = sp.simplify(alpha_derived - alpha_claimed) == 0
    print(f"[{'PASS' if ok4 else 'FAIL'}] alpha(J,m,r0) matches eq. (ads3-alpha)")

    # neutral limit J->0 and boundary limit r0->infty
    C2_holo = sp.Rational(1, 2) * sp.Abs(alpha_claimed)  # lambda=1 units
    C2_neutral = sp.simplify(C2_holo.subs(J, 0))
    C2_neutral_claimed = sp.Rational(1, 2) * m * r0 / sp.sqrt(1 + r0**2)
    ok5 = sp.simplify(C2_neutral - C2_neutral_claimed) == 0
    print(f"[{'PASS' if ok5 else 'FAIL'}] J->0 limit of C_2^holo/lambda matches m r0/(2 sqrt(1+r0^2))")

    lim = sp.limit(C2_neutral_claimed, r0, sp.oo)
    ok6 = sp.simplify(lim - m / 2) == 0
    print(f"[{'PASS' if ok6 else 'FAIL'}] r0->infty boundary limit equals m/2 (units lambda=1)")

    all_ok = all([ok1, ok2, ok3, ok4, ok5, ok6])
    print(f"\nAll symbolic checks passed: {all_ok}\n")
    return all_ok, (r, rdot, m, J, R, Pr, Energy)


# ----------------------------------------------------------------------
# 2. Numerical geodesic integration and short-time slope extraction
# ----------------------------------------------------------------------

def build_numeric_functions(r, rdot, m, J, R, Pr):
    dPr_dr = sp.diff(Pr, r)
    dPr_drdot = sp.diff(Pr, rdot)
    dR_dr = sp.diff(R, r)
    rddot_expr = sp.simplify((dR_dr - dPr_dr * rdot) / dPr_drdot)

    f_rddot = sp.lambdify((r, rdot, m, J), rddot_expr, 'numpy')
    f_Pr = sp.lambdify((r, rdot, m, J), Pr, 'numpy')
    return f_rddot, f_Pr


def alpha_analytic(m_, J_, r0_):
    return (J_**2 - m_**2 * r0_**4) / (r0_**2 * np.sqrt((m_**2 * r0_**2 + J_**2) * (1 + r0_**2)))


def naive_intercept(J_, r0_):
    return abs(J_) / r0_  # |J|/sqrt(h(r0)), h(r)=r^2


def integrate_case(f_rddot, m_, J_, r0_, t_max=5e-3):
    def rhs(t, y):
        rr, vv = y
        return [vv, f_rddot(rr, vv, m_, J_)]
    sol = solve_ivp(rhs, [0, t_max], [r0_, 0.0], method='DOP853',
                     rtol=1e-13, atol=1e-15, dense_output=True)
    return sol


def numerical_validation(r, rdot, m, J, R, Pr):
    print("=== Numerical validation (DOP853, rtol=atol~1e-13) ===")
    f_rddot, f_Pr = build_numeric_functions(r, rdot, m, J, R, Pr)

    test_cases = [
        (1.0, 0.5, 2.0),
        (1.0, 2.0, 1.0),
        (1.0, 0.1, 3.0),
        (2.0, 5.0, 0.8),
        (1.0, 0.0, 1.5),
    ]

    print(f"{'m':>5} {'J':>5} {'r0':>5} {'alpha_analytic':>16} {'alpha_numeric':>16} {'rel_err':>12}")
    rel_errs = []
    for m_, J_, r0_ in test_cases:
        sol = integrate_case(f_rddot, m_, J_, r0_)
        tw = np.linspace(1e-6, 2e-4, 200)
        ys = sol.sol(tw)
        Pw = f_Pr(ys[0], ys[1], m_, J_)
        alpha_num = np.sum(tw * Pw) / np.sum(tw * tw)
        alpha_an = alpha_analytic(m_, J_, r0_)
        relerr = abs(alpha_num - alpha_an) / abs(alpha_an) if alpha_an != 0 else abs(alpha_num)
        rel_errs.append(relerr)
        print(f"{m_:5.2f} {J_:5.2f} {r0_:5.2f} {alpha_an:16.10f} {alpha_num:16.10f} {relerr:12.3e}")

    worst = max(rel_errs)
    print(f"\nWorst relative error over the tested cases: {worst:.2e}\n")
    return f_rddot, f_Pr, worst


# ----------------------------------------------------------------------
# 3. Reproduce the validation figure referenced in the report
# ----------------------------------------------------------------------

def make_figure(f_rddot, f_Pr, out="holographic_krylov_validation.png"):
    m_ = 1.0
    r0_ = 1.5
    J_values = [0.0, 1.0, 2.0, 5.0]

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.2))

    for J_ in J_values:
        sol = integrate_case(f_rddot, m_, J_, r0_, t_max=3.6)
        t_plot = np.linspace(1e-6, 3.6, 2000)
        ys = sol.sol(t_plot)
        r_t, rdot_t = ys
        P_R = f_Pr(r_t, rdot_t, m_, J_)
        P_naive = np.sqrt((P_R / (1 / (1 + r_t**2)))**2 + (J_**2 / r_t**2)) * np.sign(P_R + 1e-30)
        # naive full-momentum magnitude |P_full| = sqrt(P_r^2/g + J^2/h) with the *original* P_r (unreduced);
        # for display purposes we track its magnitude, not the reduced P_r^(R).
        axes[0].plot(t_plot, P_R, label=f"J={J_:g}: Routhian")
        axes[0].plot(t_plot, np.full_like(t_plot, naive_intercept(J_, r0_)) if J_ != 0 else np.zeros_like(t_plot),
                     '--', color=axes[0].lines[-1].get_color(), alpha=0.5,
                     label=f"J={J_:g}: naive intercept" if J_ != 0 else None)

    axes[0].set_xlabel("t")
    axes[0].set_ylabel("momentum magnitude")
    axes[0].set_title("Naive intercept vs Routhian $P_r^{(R)}(t)$")
    axes[0].axhline(0, color='k', lw=0.5)
    axes[0].legend(fontsize=7, ncol=2)

    tw = np.linspace(2e-4, 0.05, 400)
    for J_ in J_values:
        sol = integrate_case(f_rddot, m_, J_, r0_, t_max=0.06)
        ys = sol.sol(tw)
        P_R = f_Pr(ys[0], ys[1], m_, J_)
        axes[1].plot(tw, P_R / tw, label=f"J={J_:g}")
        axes[1].axhline(alpha_analytic(m_, J_, r0_), ls=':', color=axes[1].lines[-1].get_color())
    axes[1].set_xlabel("t")
    axes[1].set_ylabel(r"$P_r^{(R)}(t)/t$")
    axes[1].set_title(r"Convergence to $\alpha$ (dotted: analytic)")
    axes[1].legend(fontsize=8)

    for J_ in J_values:
        sol = integrate_case(f_rddot, m_, J_, r0_, t_max=0.06)
        ys = sol.sol(tw)
        P_R = f_Pr(ys[0], ys[1], m_, J_)
        C_over_t2 = np.abs(P_R) / (2 * tw)  # C(t)/t^2 = (lambda/2)|alpha| for lambda=1, dot C=lambda|P_R|
        axes[2].plot(tw, C_over_t2, label=f"J={J_:g}")
        axes[2].axhline(abs(alpha_analytic(m_, J_, r0_)) / 2, ls=':', color=axes[2].lines[-1].get_color())
    axes[2].set_xlabel("t")
    axes[2].set_ylabel(r"$C_K(t)/t^2$ ($\lambda=1$)")
    axes[2].set_title(r"Convergence to $C_2$ (dotted: analytic)")
    axes[2].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out, dpi=160)
    print(f"Figure written to {out}")


if __name__ == "__main__":
    ok, (r, rdot, m, J, R, Pr, Energy) = symbolic_checks()
    f_rddot, f_Pr, worst_relerr = numerical_validation(r, rdot, m, J, R, Pr)
    make_figure(f_rddot, f_Pr)
    print(f"SUMMARY: symbolic checks passed = {ok}, worst numeric relative error = {worst_relerr:.2e}")
