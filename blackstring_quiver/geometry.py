"""
Background geometry for the N=(0,4) "black string chain" AdS3 x S^2 x CY2 x I
solutions of Couzens, Lozano, Petri & Vandoren, arXiv:2109.10413 ("N=(0,4)
Black String Chains"), specialised to the p=1 fractional-brane sector of
Section 3 of that paper.

This is the QUADRATIC generalisation, requested by the supervisors, of the
linear-rank-function quivers of Fatemiabhari, Nastase, Nunez & Roychowdhury,
arXiv:2512.14812 ("Holographic Krylov Complexity for Conformal Quiver Gauge
Theories"): there h4(eta), h8(eta) are piecewise LINEAR rank functions
(their eqs. 2.4-2.5); here, once fractional D4-branes (sourced by a
Calabi-Yau two-form H2 = gamma*omega) are switched on, h4(z) becomes
piecewise QUADRATIC while h8 = gamma is forced constant (the p=1 case of
arXiv:2109.10413, eq. (3.14)).

TAKEN FROM THE REFERENCE (arXiv:2109.10413):
  - Near-horizon Einstein-frame metric, eq. (2.11):
        ds_E^2 = e^{-Phi/2} [ (h4 h8)^{-1/2} ( ds^2_AdS3 + (1/4) ds^2_S2 )
                              + (h4/h8)^{1/2} ds^2_CY2
                              + (h4 h8)^{1/2} dz^2 ]
        e^{-Phi} = h4^{1/4} h8^{5/4}                                  (2.12)
    (the supervisors instructed us to drop the CY2 factor entirely: the
    particle does not move on it, and it never enters f1,f2,f3 below).
  - h4 is piecewise quadratic (eq. 3.31), on segments [2*pi*k, 2*pi*(k+1)],
    k = 0,...,P, with local coordinate x = z - 2*pi*k:
        h4^[k](x) = (2*pi)^2 alpha_k + 2*pi*beta_k*x - (gamma/2) x^2
    fixed by: alpha_0 = 0 (root at z=0); the continuity/recursion
        alpha_{k+1} = alpha_k + beta_k - gamma/2                       (3.32)-(3.33)
    and the closing condition (root at z=2*pi*(P+1)):
        alpha_P + beta_P - gamma/2 = 0  <=>  sum_{i=0}^{P} (beta_i - gamma/2) = 0   (3.35)-(3.36)
    Flux/SUSY quantisation forces alpha_k, beta_k, gamma all integers, gamma
    even, beta_0 >= beta_1 >= ... >= beta_P (no anti-D4/D4 mixing), and
    alpha_k > 0 for 1<=k<=P (h4 strictly positive in the interior). This
    paper's own worked example (Fig. 1, P=6) is
        beta = (22, 21, 19, 17, 15, 15), gamma = 32,
    with beta_P=beta_6 FIXED by the closing condition, not freely chosen.
  - The p=1 flux-quantisation condition gamma = h8 (eq. 3.13-3.14): we take
    this simplest sector, so h8 is the SAME integer as gamma, constant
    along the whole interval.

DERIVED HERE (straightforward algebra from the metric above, not spelled
out explicitly in the reference): restricting to a geodesic with theta =
pi/2 fixed on the S^2 (so dOmega_2^2 -> dphi^2) and dropping the CY2 and
AdS3-transverse (x_1) directions exactly as in Sec. 6 of Notes.tex, the
induced worldline Lagrangian takes the SAME form
    L = -m sqrt( F1(z)(e^{-lambda r}-rdot^2) - F2(z) zdot^2 - F3(z) phidot^2 )
as eq. (6dLagrangian) of Notes.tex, with (using e^{-Phi/2}=h4^{1/8}h8^{5/8})
    F1(z) = e^{-Phi/2}/sqrt(h4 h8)     = h4^{-3/8} h8^{1/8}
    F3(z) = F1(z)/4                     (the explicit 1/4 in front of dOmega_2^2)
    F2(z) = e^{-Phi/2} sqrt(h4 h8)     = h4^{5/8}  h8^{9/8}
so that, with h8=gamma constant,
    F1(z) = gamma^{1/8} h4(z)^{-3/8},   F3(z) = F1(z)/4,   F2(z) = gamma^{9/8} h4(z)^{5/8}.
A notable consequence, used in routhian_check.py, is that F1(z)/F3(z) = 4
identically (independent of z, unlike Quiver 1/Quiver 2 of the 6d SCFT
case): this means the J=0 equilibrium eta_max = argmax(F1) = argmax(h4) is
EXACTLY independent of J for this background (see the derivation in
Notes.tex Sec. "Extension to the quadratic black-string-chain quiver").
"""

import numpy as np

EPS_REG = 1e-12


def build_alpha(P, gamma, beta_free):
    """Given P, gamma (even integer) and the P free betas beta_0..beta_{P-1},
    returns the full beta list (length P+1, with beta_P fixed by the closing
    condition) and the full alpha list (length P+1, alpha_0=0).

    Closing condition (eq. 3.36): sum_{i=0}^{P} (beta_i - gamma/2) = 0.
    Recursion (eq. 3.32/3.33): alpha_{k+1} = alpha_k + beta_k - gamma/2.
    """
    beta_free = list(beta_free)
    if len(beta_free) != P:
        raise ValueError(f"need exactly P={P} free betas (beta_0..beta_{{P-1}}); got {len(beta_free)}")
    beta_P = gamma / 2.0 * (P + 1) - sum(beta_free)
    beta = beta_free + [beta_P]

    alpha = [0.0]
    for k in range(P):
        alpha.append(alpha[k] + beta[k] - gamma / 2.0)

    closure_residual = alpha[P] + beta[P] - gamma / 2.0
    if abs(closure_residual) > 1e-8:
        raise RuntimeError(f"closing condition not satisfied: residual={closure_residual}")

    return beta, alpha


def validate_quiver(P, gamma, beta, alpha, label=""):
    """Checks the positivity / monotonicity constraints of Sec. 3.2 of
    arXiv:2109.10413 and reports (does not raise) on any violation, since
    some of them (e.g. strict beta monotonicity) are physical-solution
    requirements rather than requirements for h4 to be a well-defined
    (if unphysical) function."""
    msgs = []
    for k in range(1, len(beta)):
        if beta[k] > beta[k - 1] + 1e-9:
            msgs.append(f"beta_{k}={beta[k]} > beta_{k-1}={beta[k-1]} (violates beta monotonicity, "
                         f"would require an anti-D4 brane)")
    for k in range(1, P + 1):
        if alpha[k] <= 0:
            msgs.append(f"alpha_{k}={alpha[k]} <= 0 (h4 would be non-positive at node {k})")
    if gamma % 2 != 0:
        msgs.append(f"gamma={gamma} is not an even integer (required by (3.33)+integrality)")
    if msgs:
        print(f"[validate_quiver{(' ' + label) if label else ''}] WARNINGS:")
        for m in msgs:
            print("   -", m)
    else:
        print(f"[validate_quiver{(' ' + label) if label else ''}] all constraints satisfied "
              f"(beta monotonic non-increasing, alpha_k>0 interior, gamma even).")
    return len(msgs) == 0


class BlackStringQuiver:
    """A single choice of (P, gamma, beta) defining h4(z), h8=gamma on
    z in [0, 2*pi*(P+1)]."""

    def __init__(self, P, gamma, beta_free, label="quiver"):
        self.P = P
        self.gamma = float(gamma)
        self.beta, self.alpha = build_alpha(P, gamma, beta_free)
        self.label = label
        self.z_max = 2 * np.pi * (P + 1)
        self.kinks = np.array([2 * np.pi * k for k in range(P + 2)])  # incl. endpoints
        validate_quiver(P, self.gamma, self.beta, self.alpha, label=label)

    def segment_index(self, z):
        k = np.floor(z / (2 * np.pi)).astype(int)
        return np.clip(k, 0, self.P)

    def h4(self, z):
        z = np.atleast_1d(np.asarray(z, dtype=float))
        k = self.segment_index(z)
        x = z - 2 * np.pi * k
        a = np.array(self.alpha)[k]
        b = np.array(self.beta)[k]
        return (2 * np.pi) ** 2 * a + 2 * np.pi * b * x - (self.gamma / 2.0) * x ** 2

    def h8(self, z):
        z = np.atleast_1d(np.asarray(z, dtype=float))
        return np.full_like(z, self.gamma)

    def metric_functions(self):
        """Returns (F1, F2, F3, e^{-Phi/2}) callables of z, matching the
        interface expected by quiver_comparison/routhian.py's f1,f2,f3."""
        gamma = self.gamma

        def h4_safe(z):
            return np.maximum(self.h4(z), EPS_REG)

        def dilaton_half(z):
            # e^{-Phi/2} = h4^{1/8} h8^{5/8}, h8=gamma constant
            return h4_safe(z) ** 0.125 * gamma ** 0.625

        def F1(z):
            return gamma ** 0.125 * h4_safe(z) ** (-0.375)

        def F3(z):
            return F1(z) / 4.0

        def F2(z):
            return gamma ** 1.125 * h4_safe(z) ** 0.625

        return F1, F2, F3, dilaton_half


# ----------------------------------------------------------------------
# The reference's own worked example (Fig. 1 of arXiv:2109.10413):
# P=6, beta_0..beta_5 = (22,21,19,17,15,15), gamma=32 (beta_6=3 fixed by closure).
# TAKEN FROM THE REFERENCE.
# ----------------------------------------------------------------------
REFERENCE_EXAMPLE = dict(P=6, gamma=32, beta_free=[22, 21, 19, 17, 15, 15])


def make_reference_quiver():
    return BlackStringQuiver(**REFERENCE_EXAMPLE, label="arXiv:2109.10413 Fig.1 example")


def _validate():
    print("=== Black-string-chain geometry self-validation ===")
    q = make_reference_quiver()
    print("beta:", q.beta)
    print("alpha:", q.alpha)

    z = np.linspace(1e-3, q.z_max - 1e-3, 20001)
    h4v = q.h4(z)
    print(f"h4 range: [{h4v.min():.4f}, {h4v.max():.4f}] (min should be > 0 strictly in interior)")
    print(f"h4(0)={q.h4(np.array([1e-8]))[0]:.3e}, h4(z_max)={q.h4(np.array([q.z_max-1e-8]))[0]:.3e} "
          f"(both should be ~0)")

    # continuity check at each interior kink (h4 continuous but h4' need not be)
    for k in range(1, q.P + 1):
        bp = 2 * np.pi * k
        left = q.h4(np.array([bp - 1e-7]))[0]
        right = q.h4(np.array([bp + 1e-7]))[0]
        print(f"  continuity at kink k={k} (z={bp:.4f}): h4(bp-)={left:.6f}, h4(bp+)={right:.6f}, "
              f"diff={abs(left-right):.2e}")

    F1, F2, F3, edPhi = q.metric_functions()
    print(f"F1/F3 identically 4 check: max|F1/F3-4| = {np.max(np.abs(F1(z)/F3(z)-4)):.3e} (must be 0)")

    from scipy.signal import argrelextrema
    maxima = argrelextrema(h4v, np.greater)[0]
    minima = argrelextrema(h4v, np.less)[0]
    print(f"h4 has {len(maxima)} local maxima (interior 'wells' of the F1-potential) "
          f"at z = {np.round(z[maxima], 3)}")
    print(f"h4 has {len(minima)} local minima (interior 'barriers') at z = {np.round(z[minima], 3)} "
          f"(should coincide with kinks 2*pi*k for k=1..P)")
    print(f"global max of h4 (the J-independent equilibrium of the eta-motion) at "
          f"z = {z[np.argmax(h4v)]:.4f}, in segment k={q.segment_index(np.array([z[np.argmax(h4v)]]))[0]}")


if __name__ == "__main__":
    _validate()
