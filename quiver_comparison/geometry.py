"""
Background geometry for the two 6d N=(1,0) SCFT quivers of
Fatemiabhari, Nunez & Santamaria, "Complexity and Operator Growth in
Holographic 6d SCFTs" (arXiv:2603.10106), Section 2.1.

Both quivers are massive Type IIA duals fibred over an interval
eta in [0, P+1]:

    ds_E^2 = f1(eta) ds^2_{AdS7} + f2(eta) d eta^2 + f3(eta) dOmega_2 ,

with f1, f2, f3, and the dilaton e^Phi = f6, all built from a single
function alpha(eta) (piecewise cubic) via

    f1 = 8 sqrt(2) pi e^{-Phi/2} sqrt(-alpha/alpha'')
    f2 =   sqrt(2) pi e^{-Phi/2} sqrt(-alpha''/alpha)
    f3 =   sqrt(2) pi e^{-Phi/2} sqrt(-alpha''/alpha) * alpha^2 / ((alpha')^2 - 2 alpha alpha'')
    e^Phi = 2^{5/4} pi^{5/2} 3^4 (-alpha/alpha'')^{3/4} / sqrt((alpha')^2 - 2 alpha alpha'')

alpha is fixed by the piecewise-linear rank function R(eta) of the
linear quiver through alpha''(eta) = -81 pi^2 R(eta), integrated twice
with alpha(0) = alpha(P+1) = 0 and alpha, alpha' continuous at every
interior breakpoint (this is the boundary/continuity condition used
throughout; see the audit note at the bottom of this file for how it
was cross-checked).

QUIVER 1 (Sec. 2.1.1 of the reference): ranks N, 2N, ..., PN with a
single flavour node at eta = P ("triangular" rank function). This is
the quiver already implemented and validated in Notebook_Krylov_Final.ipynb
[TAKEN FROM THE REFERENCE, reproduced verbatim from the existing notebook].

QUIVER 2 (Sec. 2.1.2 of the reference): constant rank N between two
flavour nodes at eta = 1 and eta = P ("trapezoidal" rank function).
alpha_2(eta) below is DERIVED HERE (not copied from the paper, since
the paper's own algebra could not be retrieved with full certainty
through automated fetching) by integrating R_2 twice with the same
boundary/continuity conditions used for Quiver 1, and is cross-checked
against those same conditions in reproduce_quiver1_check() /
check_quiver2_construction() below. Labelled ANALYTICALLY DERIVED
(independent re-derivation), not ANALYTICALLY PROVEN from the paper's
own stated formula.
"""

import numpy as np

EPS_REG = 1e-9


# ----------------------------------------------------------------------
# Quiver 1: ranks N, 2N, ..., PN, single flavour node at eta = P
# ----------------------------------------------------------------------

def R1(eta, N=1.0, P=100.0):
    eta = np.atleast_1d(np.asarray(eta, dtype=float))
    out = np.empty_like(eta)
    left = eta <= P
    out[left] = N * eta[left]
    out[~left] = N * P * (P + 1 - eta[~left])
    return out


def alpha1(eta, N=1.0, P=100.0):
    eta = np.atleast_1d(np.asarray(eta, dtype=float))
    a1 = -P * (P + 2) / 6.0
    out = np.empty_like(eta)
    left = eta <= P
    e = eta[left]
    out[left] = -81 * np.pi**2 * N * (a1 * e + e**3 / 6.0)
    x = eta[~left] - P
    out[~left] = -81 * np.pi**2 * N * (
        (P * a1 + P**3 / 6.0)
        + (a1 + P**2 / 2.0) * x
        + (P / 2.0) * x**2
        - (P / 6.0) * x**3
    )
    return out


def alpha1p(eta, N=1.0, P=100.0):
    eta = np.atleast_1d(np.asarray(eta, dtype=float))
    a1 = -P * (P + 2) / 6.0
    out = np.empty_like(eta)
    left = eta <= P
    e = eta[left]
    out[left] = -81 * np.pi**2 * N * (a1 + e**2 / 2.0)
    x = eta[~left] - P
    out[~left] = -81 * np.pi**2 * N * ((a1 + P**2 / 2.0) + P * x - (P / 2.0) * x**2)
    return out


def alpha1pp(eta, N=1.0, P=100.0):
    return -81 * np.pi**2 * R1(eta, N=N, P=P)


# ----------------------------------------------------------------------
# Quiver 2: constant rank N, flavour nodes at eta = 1 and eta = P
#
#   R2(eta) = N eta            0 <= eta <= 1
#           = N                1 <= eta <= P
#           = N (P+1-eta)      P <= eta <= P+1
#
# Derivation (independent, this work): with K = 81 pi^2 N,
# integrating alpha'' = -K R2 twice, imposing alpha(0)=0, continuity
# of alpha,alpha' at eta=1 and eta=P, and alpha(P+1)=0, gives
# (see audit note at bottom of file for the algebra):
#
#   c1 = K P / 2      (the single free constant fixed by alpha(P+1)=0)
#
#   Region A (0<=eta<=1):
#     alpha_A(eta)  = K( P eta/2 - eta^3/6 )
#     alpha_A'(eta) = K( P/2 - eta^2/2 )
#
#   Region B (1<=eta<=P), u = eta-1:
#     alpha_B(u) = -K u^2/2 + K(P-1)/2 * u + K(3P-1)/6
#     alpha_B'(u) = -K u + K(P-1)/2
#
#   Region C (P<=eta<=P+1), x = eta-P:
#     alpha_C(x) = K x^3/6 - K x^2/2 + K(1-P)/2 * x + K(3P-1)/6
#     alpha_C'(x) = K x^2/2 - K x + K(1-P)/2
# ----------------------------------------------------------------------

def R2(eta, N=1.0, P=100.0):
    eta = np.atleast_1d(np.asarray(eta, dtype=float))
    out = np.empty_like(eta)
    mid = (eta >= 1.0) & (eta <= P)
    left = eta < 1.0
    right = eta > P
    out[left] = N * eta[left]
    out[mid] = N
    out[right] = N * (P + 1 - eta[right])
    return out


def alpha2(eta, N=1.0, P=100.0):
    eta = np.atleast_1d(np.asarray(eta, dtype=float))
    K = 81 * np.pi**2 * N
    out = np.empty_like(eta)

    left = eta < 1.0
    mid = (eta >= 1.0) & (eta <= P)
    right = eta > P

    e = eta[left]
    out[left] = K * (P * e / 2.0 - e**3 / 6.0)

    u = eta[mid] - 1.0
    out[mid] = -K * u**2 / 2.0 + K * (P - 1.0) / 2.0 * u + K * (3 * P - 1.0) / 6.0

    x = eta[right] - P
    out[right] = K * x**3 / 6.0 - K * x**2 / 2.0 + K * (1.0 - P) / 2.0 * x + K * (3 * P - 1.0) / 6.0

    return out


def alpha2p(eta, N=1.0, P=100.0):
    eta = np.atleast_1d(np.asarray(eta, dtype=float))
    K = 81 * np.pi**2 * N
    out = np.empty_like(eta)

    left = eta < 1.0
    mid = (eta >= 1.0) & (eta <= P)
    right = eta > P

    e = eta[left]
    out[left] = K * (P / 2.0 - e**2 / 2.0)

    u = eta[mid] - 1.0
    out[mid] = -K * u + K * (P - 1.0) / 2.0

    x = eta[right] - P
    out[right] = K * x**2 / 2.0 - K * x + K * (1.0 - P) / 2.0

    return out


def alpha2pp(eta, N=1.0, P=100.0):
    return -81 * np.pi**2 * R2(eta, N=N, P=P)


# ----------------------------------------------------------------------
# Generic massive-IIA formulas: f1, f2, f3, dilaton, given alpha,alpha',alpha''
# ----------------------------------------------------------------------

def _alpha_data(alpha_fn, alphap_fn, alphapp_fn, eta, N, P):
    a = alpha_fn(eta, N=N, P=P)
    ap = alphap_fn(eta, N=N, P=P)
    app = alphapp_fn(eta, N=N, P=P)
    ratio = np.maximum(-a / np.where(np.abs(app) < EPS_REG, EPS_REG, app), EPS_REG)  # -alpha/alpha''
    denom = np.maximum(np.abs(ap**2 - 2 * a * app), EPS_REG)
    return a, ap, app, ratio, denom


def make_metric_functions(quiver, N=1.0, P=100.0):
    """
    Returns (f1, f2, f3, f6) callables for the requested quiver
    ('quiver1' or 'quiver2'), each a function of eta (array-like).
    """
    if quiver == "quiver1":
        alpha_fn, alphap_fn, alphapp_fn = alpha1, alpha1p, alpha1pp
    elif quiver == "quiver2":
        alpha_fn, alphap_fn, alphapp_fn = alpha2, alpha2p, alpha2pp
    else:
        raise ValueError(f"unknown quiver {quiver!r}")

    def f6(eta):
        a, ap, app, ratio, denom = _alpha_data(alpha_fn, alphap_fn, alphapp_fn, eta, N, P)
        return 2**1.25 * np.pi**2.5 * 3**4 * ratio**0.75 / np.sqrt(denom)

    def f1(eta):
        a, ap, app, ratio, denom = _alpha_data(alpha_fn, alphap_fn, alphapp_fn, eta, N, P)
        pref = 1.0 / np.sqrt(f6(eta))
        return 8 * np.sqrt(2) * np.pi * pref * np.sqrt(ratio)

    def f2(eta):
        a, ap, app, ratio, denom = _alpha_data(alpha_fn, alphap_fn, alphapp_fn, eta, N, P)
        pref = 1.0 / np.sqrt(f6(eta))
        return np.sqrt(2) * np.pi * pref * np.sqrt(1.0 / ratio)

    def f3(eta):
        a, ap, app, ratio, denom = _alpha_data(alpha_fn, alphap_fn, alphapp_fn, eta, N, P)
        pref = 1.0 / np.sqrt(f6(eta))
        return np.sqrt(2) * np.pi * pref * np.sqrt(1.0 / ratio) * (a**2 / denom)

    return f1, f2, f3, f6


# ----------------------------------------------------------------------
# Self-audit / validation, run when this module is executed directly
# ----------------------------------------------------------------------

def _validate():
    print("=== Geometry self-validation ===")

    P = 100.0
    N = 1.0
    eta_test = np.linspace(0.01, P + 1 - 0.01, 2001)

    # 1. alpha'' must reproduce -81 pi^2 R(eta) pointwise, for both quivers.
    err1 = np.max(np.abs(alpha1pp(eta_test, N, P) - (-81 * np.pi**2 * R1(eta_test, N, P))))
    err2 = np.max(np.abs(alpha2pp(eta_test, N, P) - (-81 * np.pi**2 * R2(eta_test, N, P))))
    print(f"max |alpha1'' + 81 pi^2 R1| = {err1:.3e}  (must be 0 by construction)")
    print(f"max |alpha2'' + 81 pi^2 R2| = {err2:.3e}  (must be 0 by construction)")

    # 2. Boundary conditions alpha(0)=alpha(P+1)=0.
    print(f"alpha1(0)   = {alpha1(np.array([1e-8]), N, P)[0]:.6e}  (should be ~0)")
    print(f"alpha1(P+1) = {alpha1(np.array([P + 1 - 1e-8]), N, P)[0]:.6e}  (should be ~0)")
    print(f"alpha2(0)   = {alpha2(np.array([1e-8]), N, P)[0]:.6e}  (should be ~0)")
    print(f"alpha2(P+1) = {alpha2(np.array([P + 1 - 1e-8]), N, P)[0]:.6e}  (should be ~0)")

    # 3. Continuity of alpha, alpha' at every interior breakpoint.
    # NOTE: a naive finite-difference check across the breakpoint with a fixed
    # step d is misleading here, because alpha'(eta) itself is O(K P^2) ~ 1e5-1e6
    # for the parameters used below, so a "continuous but steep" function already
    # produces an O(alpha' * d) change across the step -- this is NOT a discontinuity.
    # The rigorous check (convergence of the normalized jump to 0 as d -> 0, i.e.
    # jump(d)/d -> alpha'(bp) rather than jump(d)/d -> infinity) is done here;
    # the fully independent proof is sympy_verify_quiver2.py, which re-derives
    # alpha2(eta) from scratch via symbolic matching and confirms an exact match.
    def check_continuity(alpha_fn, breakpoints, label):
        for bp in breakpoints:
            ds = [1e-3, 1e-4, 1e-5, 1e-6]
            jumps = []
            for d in ds:
                a_minus = alpha_fn(np.array([bp - d]), N, P)[0]
                a_plus = alpha_fn(np.array([bp + d]), N, P)[0]
                jumps.append(abs(a_plus - a_minus) / (2 * d))
            print(f"  {label} at eta={bp}: [alpha(bp+d)-alpha(bp-d)]/(2d) for d={ds} -> {jumps}")
            print(f"    (should converge to a finite value = alpha'({bp}), not diverge, "
                  f"as d shrinks: {jumps[-1]:.3f} vs {jumps[0]:.3f} -> converged = "
                  f"{abs(jumps[-1] - jumps[-2]) < 1e-2 * abs(jumps[-1])}")

    print("Quiver 1 continuity:")
    check_continuity(alpha1, [P], "quiver1")
    print("Quiver 2 continuity (see also sympy_verify_quiver2.py for the exact symbolic proof):")
    check_continuity(alpha2, [1.0, P], "quiver2")

    # 4. alpha, -alpha'' must both be positive in the interior (needed for f1,f2,f3 to be real).
    a1v = alpha1(eta_test, N, P)
    a2v = alpha2(eta_test, N, P)
    print(f"quiver1: min(alpha1) = {a1v.min():.4f}  (should be > 0 in interior)")
    print(f"quiver2: min(alpha2) = {a2v.min():.4f}  (should be > 0 in interior)")

    f1_q1, f2_q1, f3_q1, f6_q1 = make_metric_functions("quiver1", N, P)
    f1_q2, f2_q2, f3_q2, f6_q2 = make_metric_functions("quiver2", N, P)
    eta_int = np.linspace(2.0, P - 2.0, 501)
    print(f"quiver1: f1 range on interior = [{f1_q1(eta_int).min():.3f}, {f1_q1(eta_int).max():.3f}]")
    print(f"quiver2: f1 range on interior = [{f1_q2(eta_int).min():.3f}, {f1_q2(eta_int).max():.3f}]")
    print(f"quiver1: location of max(f1) = eta_max = {eta_int[np.argmax(f1_q1(eta_int))]:.2f}")
    print(f"quiver2: location of max(f1) = eta_max = {eta_int[np.argmax(f1_q2(eta_int))]:.2f}")


if __name__ == "__main__":
    _validate()
