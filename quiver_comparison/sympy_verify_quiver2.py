"""
Independent symbolic (SymPy) re-derivation of alpha_2(eta) for Quiver 2,
as a cross-check on the hand derivation used in geometry.py.

Procedure: alpha'' = -81 pi^2 R2(eta) on each of the three pieces
[0,1], [1,P], [P,P+1]; integrate twice per piece (2 constants each,
6 unknowns total); impose alpha(0)=0, alpha(P+1)=0, and continuity of
alpha and alpha' at eta=1 and eta=P (6 equations); solve; compare
the result to the closed-form expressions implemented in geometry.py.
"""

import sympy as sp

eta, N, P, K = sp.symbols('eta N P K', positive=True)

# --- piece A: 0 <= eta <= 1 ---
cA1, cA0 = sp.symbols('cA1 cA0')
alphaA = sp.integrate(sp.integrate(-K * (N * eta / N), eta) + cA1, eta) + cA0
# (R2 = N*eta on this piece; write alpha'' = -K*eta after dividing out N, i.e. K = 81 pi^2 N)
alphaA = -K * eta**3 / 6 + cA1 * eta + cA0

# --- piece B: 1 <= eta <= P --- (R2 = N, alpha'' = -K)
cB1, cB0 = sp.symbols('cB1 cB0')
alphaB = -K * eta**2 / 2 + cB1 * eta + cB0

# --- piece C: P <= eta <= P+1 --- (R2 = N*(P+1-eta), alpha'' = -K*(P+1-eta))
cC1, cC0 = sp.symbols('cC1 cC0')
alphaC = K * eta**3 / 6 - K * (P + 1) * eta**2 / 2 + cC1 * eta + cC0

alphaAp = sp.diff(alphaA, eta)
alphaBp = sp.diff(alphaB, eta)
alphaCp = sp.diff(alphaC, eta)

eqs = [
    sp.Eq(alphaA.subs(eta, 0), 0),                      # alpha(0)=0
    sp.Eq(alphaA.subs(eta, 1), alphaB.subs(eta, 1)),    # continuity of alpha at eta=1
    sp.Eq(alphaAp.subs(eta, 1), alphaBp.subs(eta, 1)),  # continuity of alpha' at eta=1
    sp.Eq(alphaB.subs(eta, P), alphaC.subs(eta, P)),    # continuity of alpha at eta=P
    sp.Eq(alphaBp.subs(eta, P), alphaCp.subs(eta, P)),  # continuity of alpha' at eta=P
    sp.Eq(alphaC.subs(eta, P + 1), 0),                  # alpha(P+1)=0
]

sol = sp.solve(eqs, [cA1, cA0, cB1, cB0, cC1, cC0], dict=True)[0]
print("Solved integration constants (SymPy, independent of geometry.py):")
for k, v in sol.items():
    print(f"  {k} = {sp.simplify(v)}")

alphaA_sol = sp.simplify(alphaA.subs(sol))
alphaB_sol = sp.simplify(alphaB.subs(sol))
alphaC_sol = sp.simplify(alphaC.subs(sol))

print("\nalpha_A(eta) [0<=eta<=1]   =", alphaA_sol)
print("alpha_B(eta) [1<=eta<=P]   =", sp.expand(alphaB_sol))
print("alpha_C(eta) [P<=eta<=P+1] =", sp.expand(alphaC_sol))

# --- Compare against geometry.py's closed form ---
# geometry.py:
#   region A: K*(P*eta/2 - eta**3/6)
#   region B (u = eta-1): -K*u**2/2 + K*(P-1)/2*u + K*(3P-1)/6
#   region C (x = eta-P): K*x**3/6 - K*x**2/2 + K*(1-P)/2*x + K*(3P-1)/6

u = eta - 1
x = eta - P
alphaA_geom = K * (P * eta / 2 - eta**3 / 6)
alphaB_geom = -K * u**2 / 2 + K * (P - 1) / 2 * u + K * (3 * P - 1) / 6
alphaC_geom = K * x**3 / 6 - K * x**2 / 2 + K * (1 - P) / 2 * x + K * (3 * P - 1) / 6

diffA = sp.simplify(alphaA_sol - alphaA_geom)
diffB = sp.simplify(alphaB_sol - alphaB_geom)
diffC = sp.simplify(alphaC_sol - alphaC_geom)

print("\n=== Cross-check against geometry.py closed form ===")
print("alpha_A(sympy) - alpha_A(geometry.py) =", diffA, "  (must be 0)")
print("alpha_B(sympy) - alpha_B(geometry.py) =", diffB, "  (must be 0)")
print("alpha_C(sympy) - alpha_C(geometry.py) =", diffC, "  (must be 0)")

assert diffA == 0 and diffB == 0 and diffC == 0, "MISMATCH: geometry.py alpha2 does not match independent SymPy re-derivation!"
print("\nPASSED: geometry.py's alpha2(eta) matches an independent SymPy re-derivation exactly.")

# --- Also verify alpha'' reproduces -K*R2 on each piece ---
R2A, R2B, R2C = eta, sp.Integer(1), (P + 1 - eta)
print("\nalpha_A'' + K*R2_A =", sp.simplify(sp.diff(alphaA_sol, eta, 2) + K * R2A))
print("alpha_B'' + K*R2_B =", sp.simplify(sp.diff(alphaB_sol, eta, 2) + K * R2B))
print("alpha_C'' + K*R2_C =", sp.simplify(sp.diff(alphaC_sol, eta, 2) + K * R2C))
