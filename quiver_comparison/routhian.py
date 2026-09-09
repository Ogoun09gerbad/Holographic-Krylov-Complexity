"""
Fixed-charge Routhian dynamics for a massive geodesic probe in the 6d SCFT
background of Fatemiabhari, Nunez & Santamaria (arXiv:2603.10106), Sec. 3,
as derived in Section 6 of Notes.tex ("Extension: Fixed-charge Routhian for
holographic 6d SCFT backgrounds"). Everything here is expressed generically
in terms of the metric functions f1(eta), f2(eta), f3(eta) of a given
quiver, so the SAME code applies unchanged to Quiver 1 or Quiver 2 -- this
is the numerical embodiment of the claim that the Routhian mechanism is
geometry-independent.

Lagrangian (m=1 units, as in the reference and the existing notebook):
    L = sqrt( f1(eta)(e^{-2r} - rdot^2) - f2(eta) etadot^2 - f3(eta) phidot^2 )

Conserved quantities:
    J = f3(eta) phidot / D                (R-charge; D = the sqrt above)
    H = f1(eta) e^{-2r} / D                (Hamiltonian; conserved since
                                             the Lagrangian has no explicit t)

Fixed-charge Routhian (eq. 6.20 of Notes.tex):
    R(r,eta,rdot,etadot;J) = -m_eff(eta) sqrt( f1(eta)(e^{-2r}-rdot^2) - f2(eta) etadot^2 )
    m_eff(eta) = sqrt( 1 + J^2/f3(eta) )        (m=1)

Radial motion decouples EXACTLY (eq. 6.??, "An exact all-t check" in
Notes.tex): for a probe released from rest at r(0)=r_UV,
    r(t)    = 0.5 * log( e^{2 r_UV} + t^2 )
    P_r(t)  = H * t                              (exact, for ANY J and ANY eta(t))

so only the quiver-direction motion eta(t) needs to be integrated
numerically, exactly as in the existing notebook.
"""

import numpy as np
from scipy.integrate import solve_ivp


def m_eff2(eta, J, f3):
    return 1.0 + J**2 / f3(eta)


def r_of_t(t, r_UV):
    return 0.5 * np.log(np.exp(2 * r_UV) + t**2)


def rdot_of_t(t, r_UV):
    return t / (np.exp(2 * r_UV) + t**2)


def H_of_J(J, eta0, r_UV, f1, f3):
    """Conserved Hamiltonian, fixed by release-from-rest initial data
    r(0)=r_UV, eta(0)=eta0, rdot(0)=etadot(0)=0."""
    eta0_arr = np.array([eta0])
    return np.sqrt(f1(eta0_arr))[0] * np.sqrt(m_eff2(eta0_arr, J, f3))[0] * np.exp(-r_UV)


def etadot_of_Peta(eta, Peta, t, J, r_UV, f1, f2, f3):
    """Invert P_eta = f2(eta) etadot / D_R for etadot. Smooth through
    etadot=0 (no sign ambiguity, unlike solving the first-order 'energy
    conservation' equation directly for etadot)."""
    rt, rdt = r_of_t(t, r_UV), rdot_of_t(t, r_UV)
    kin = np.maximum(f1(eta) * (np.exp(-2 * rt) - rdt**2), 0.0)
    denom = f2(eta) * (m_eff2(eta, J, f3) * f2(eta) + Peta**2)
    return Peta * np.sqrt(kin / denom)


def R_routhian(eta, etadot, t, J, r_UV, f1, f2, f3):
    rt, rdt = r_of_t(t, r_UV), rdot_of_t(t, r_UV)
    DR2 = np.maximum(f1(eta) * (np.exp(-2 * rt) - rdt**2) - f2(eta) * etadot**2, 1e-300)
    return -np.sqrt(m_eff2(eta, J, f3)) * np.sqrt(DR2)


def geodesic_rhs(t, y, J, r_UV, f1, f2, f3):
    """Hamilton-type equations for phase-space state (eta, P_eta)."""
    eta, Peta = np.array([y[0]]), np.array([y[1]])
    etadot = etadot_of_Peta(eta, Peta, t, J, r_UV, f1, f2, f3)
    deta = 1e-6 * max(1.0, abs(eta[0]))
    dRdeta = (
        R_routhian(eta + deta, etadot, t, J, r_UV, f1, f2, f3)
        - R_routhian(eta - deta, etadot, t, J, r_UV, f1, f2, f3)
    ) / (2 * deta)
    return [etadot[0], dRdeta[0]]


def integrate_trajectory(J, eta0, r_UV, t_max, f1, f2, f3, max_step=0.01, atol=1e-11, rtol=1e-11):
    sol = solve_ivp(
        geodesic_rhs, [0.0, t_max], [eta0, 0.0], args=(J, r_UV, f1, f2, f3),
        dense_output=True, atol=atol, rtol=rtol, max_step=max_step,
    )
    return sol


def naive_and_routhian_momenta(t_arr, sol, J, eta0, r_UV, f1, f2, f3):
    """Given an integrated (eta, P_eta) trajectory, compute the naive
    3-direction proper momentum and the fixed-charge (Routhian) 2-direction
    proper momentum at each time in t_arr."""
    eta_t, Peta_t = sol.sol(t_arr)
    rt, rdt = r_of_t(t_arr, r_UV), rdot_of_t(t_arr, r_UV)
    etadot_t = etadot_of_Peta(eta_t, Peta_t, t_arr, J, r_UV, f1, f2, f3)

    H = H_of_J(J, eta0, r_UV, f1, f3)
    D = f1(eta_t) * np.exp(-2 * rt) / H
    phidot_t = J * f1(eta_t) * np.exp(-2 * rt) / (H * f3(eta_t))

    rho_naive_dot = np.sqrt(f1(eta_t) * rdt**2 + f2(eta_t) * etadot_t**2 + f3(eta_t) * phidot_t**2)
    rho_R_dot = np.sqrt(f1(eta_t) * rdt**2 + f2(eta_t) * etadot_t**2)

    P_naive = rho_naive_dot / D
    P_R = rho_R_dot / D
    P_r = f1(eta_t) * rdt / D

    return dict(eta=eta_t, Peta=Peta_t, etadot=etadot_t, D=D, H=H,
                P_naive=P_naive, P_R=P_R, P_r=P_r, phidot=phidot_t)


# ------------------------------------------------------------------
# Analytic short-time coefficients (eqs. 6.21-6.24 of Notes.tex)
# ------------------------------------------------------------------

def A_r_analytic(J, eta0, r_UV, f1, f3):
    eta0_arr = np.array([eta0])
    return np.sqrt(m_eff2(eta0_arr, J, f3))[0] * np.sqrt(f1(eta0_arr))[0] * np.exp(-r_UV)


def A_eta_analytic(J, eta0, r_UV, f1, f3, deta=1e-4):
    eta_pm = np.array([eta0 - deta, eta0 + deta])
    V = np.sqrt(m_eff2(eta_pm, J, f3)) * np.sqrt(f1(eta_pm))
    Vprime = (V[1] - V[0]) / (2 * deta)
    return -Vprime * np.exp(-r_UV)


def short_time_coefficients(J, eta0, r_UV, f1, f2, f3):
    """Returns dict with r2, eta2, A_r, A_eta, A_6d, naive_intercept."""
    Ar = A_r_analytic(J, eta0, r_UV, f1, f3)
    Aeta = A_eta_analytic(J, eta0, r_UV, f1, f3)
    eta0_arr = np.array([eta0])
    r2 = np.exp(-2 * r_UV)
    eta2 = (Aeta * np.sqrt(f1(eta0_arr))[0] * np.exp(-r_UV)
            / (np.sqrt(m_eff2(eta0_arr, J, f3))[0] * f2(eta0_arr)[0]))
    A6d = (Ar * r2 + Aeta * eta2) / np.sqrt(f1(eta0_arr)[0] * r2**2 + f2(eta0_arr)[0] * eta2**2)
    naive_intercept = abs(J) / np.sqrt(f3(eta0_arr))[0]
    return dict(A_r=Ar, A_eta=Aeta, r2=r2, eta2=eta2, A_6d=A6d, naive_intercept=naive_intercept)
