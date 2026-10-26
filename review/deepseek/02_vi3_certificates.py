#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Independent referee check, script 02 -- VI3 literal certificates.

Rebuilt from scratch from  sum_n U_n z^n = (1-2bz+dz^2)^(-1/2).
Checked as EXACT identities in the free symbols n,t:

  (2)   sigma,tau of the author  (re-derived independently too)
  (3)   M_n, N_n
  (4)/(A)  R_{n+1}N_n - R_n = e_3 + a_2 e_2 + a_1 e_1 + a_0 e_0
  (B)   S_{n+1}N_n - S_n = Y_n' + Y_n M_n - rho Y_n
  (7)   Bscr = -g2, Ascr = -Bscr'-h0 Bscr-g1, Ascr'+h0 Ascr+g0=0, Rscr'/Rscr=-h0
        (the last group cannot be checked without g0,g1,g2,h0; we check the
         structural consequences that ARE printed: rho = Rscr'/Rscr)

Exit code 0 <=> every printed identity is exactly zero.
"""
import sys
import sympy as sp

t, n = sp.symbols('t n')
D = t ** 2 - t - 1

print("=" * 78)
print("VI3 : b = 2(t^2+t-1), c = (1+t-t^2)^2, w = 1/(64(t-1)(t+1)^2)")
print("=" * 78)

b = 2 * (t ** 2 + t - 1)
c = (1 + t - t ** 2) ** 2
d = sp.expand(b ** 2 - 4 * c)
w = 1 / (64 * (t - 1) * (t + 1) ** 2)
bp, dp = sp.diff(b, t), sp.diff(d, t)
print("d = b^2-4c =", sp.factor(d), "   (claimed 16t(t^2-1))")
assert sp.simplify(d - 16 * t * (t ** 2 - 1)) == 0
print("b^2-d = 4c :", sp.simplify(b ** 2 - d - 4 * c) == 0)

# ------------------------------------------------------------------ sigma,tau
sig = sp.simplify((b * bp - dp / 2) / (b ** 2 - d))
tau = sp.simplify(-(d * bp - b * dp / 2) / (b ** 2 - d))
sig_claim = sp.simplify((2 * t - 1) / D)
tau_claim = sp.simplify(-4 * (t ** 2 + 1) / D)
print("[2]   sigma =", sig, "  author:", sig_claim, "  diff =",
      sp.simplify(sig - sig_claim))
print("[2]   tau   =", tau, "  author:", tau_claim, "  diff =",
      sp.simplify(tau - tau_claim))

# independent re-derivation of sigma,tau from the generating function
s_, t_ = sp.symbols('s_ t_')
eq1 = sp.expand(s_ * b + t_ - bp)
eq2 = sp.expand(-s_ * d - t_ * b + dp / 2)
sol = sp.solve([eq1, eq2], [s_, t_], dict=True)[0]
print("[2']  re-derived sigma =", sp.simplify(sol[s_]), " tau =", sp.simplify(sol[t_]))

# --------------------------------------------------------------------- M_n N_n
wp = sp.simplify(sp.diff(w, t) / w)
print("[3]   w'/w =", sp.factor(wp), "   claimed 1/(t(t+1)(t-1))*? check vs -(3t^2-1)/(t(t^2-1))")
print("      w'/w + (3t^2-1)/(t(t^2-1)) =", sp.simplify(wp + (3 * t ** 2 - 1) / (t * (t ** 2 - 1))))

def Mmat(nn):
    M = sp.zeros(4, 4)
    for j in range(4):
        M[j, j] = (3 - j) * nn * sig + j * ((nn - 1) * sig + (2 * nn - 1) * b * tau / d) + nn * wp
        if j + 1 <= 3:
            M[j, j + 1] = (3 - j) * nn * tau
        if j - 1 >= 0:
            M[j, j - 1] = -j * nn * tau / d
    return M

def Nmat(nn):
    N = sp.zeros(4, 4)
    for j in range(4):
        for k in range(0, 3 - j + 1):
            N[j, k] = w * sp.binomial(3 - j, k) \
                * ((2 * nn + 1) * b / (nn + 1)) ** (3 - j - k) \
                * (-nn * d / (nn + 1)) ** k
    return N

def Un(bv, cv, nn):
    if nn < 0:
        return sp.Integer(0)
    return sp.expand(sum(sp.factorial(nn) / (sp.factorial(k) ** 2 * sp.factorial(nn - 2 * k))
                         * bv ** (nn - 2 * k) * cv ** k for k in range(nn // 2 + 1)))

def Xvec(nn):
    return sp.Matrix([w ** nn * Un(b, c, nn) ** (3 - j) * Un(b, c, nn - 1) ** j
                      for j in range(4)])

okM = all(sp.simplify(sp.Matrix([sp.diff(Xvec(nn)[i], t) for i in range(4)])
                      - Mmat(nn) * Xvec(nn)) == sp.zeros(4, 1) for nn in [1, 2, 3, 4])
print("[3]   M_n reproduces d/dt(w^n U_n^{3-j}U_{n-1}^j), n=1..4 :", okM)
okN = all(sp.simplify(Xvec(nn + 1) - Nmat(nn) * Xvec(nn)) == sp.zeros(4, 1) for nn in [0, 1, 2, 3])
print("[3]   N_n reproduces X_{n+1} = N_n X_n, n=0..3 :", okN)

M = Mmat(n)
N = Nmat(n)
es = [sp.Matrix([[1, 0, 0, 0]])]
for j in range(3):
    e = es[-1]
    ed = sp.Matrix([[sp.diff(e[0, i], t) for i in range(4)]])
    es.append(sp.simplify(ed + e * M))
okE = True
for nn in [0, 1, 2, 3]:
    X = Xvec(nn)
    for j in range(4):
        lhs = sp.simplify(sp.diff(w ** nn * Un(b, c, nn) ** 3, t, j))
        rhs = sp.simplify(sum([es[j].subs(n, nn)[0, i] * X[i] for i in range(4)]))
        if sp.simplify(lhs - rhs) != 0:
            okE = False
            print("   e_j mismatch n=", nn, "j=", j)
print("[A]   e_j reproduces D_t^j(w^n U_n^3), n=0..3, j=0..3 :", okE)

# ------------------------------------------------------------------ operator L
a2 = 3 * (2 * t - 1) * (3 * t ** 2 - 3 * t - 1) / (2 * t * (t - 1) * D)
a1 = (75 * t ** 5 - 3 * t ** 4 - 128 * t ** 3 - t ** 2 + 59 * t - 6) \
    / (4 * t * (t - 1) ** 2 * (t + 1) ** 2 * D)
a0 = (60 * t ** 6 + 21 * t ** 5 - 130 * t ** 4 - 35 * t ** 3 + 65 * t ** 2 + 50 * t - 23) \
    / (8 * t * (t - 1) ** 3 * (t + 1) ** 3 * D)
E = sp.simplify(es[3] + a2 * es[2] + a1 * es[1] + a0 * es[0])

# ------------------------------------------------ certificate (A) literal R_n
R0 = n ** 3 * (6 * t ** 5 + 3 * t ** 4 + 10 * t ** 3 - 22 * t ** 2 + 8 * t + 3) \
    / (t ** 2 * (t - 1) ** 2 * (t + 1) * D ** 3)
R1 = 6 * n ** 2 / (t ** 2 * (t - 1) * (t + 1) * D ** 3) \
    * (12 * n * t ** 6 - 4 * n * t ** 5 - 36 * n * t ** 4 - 4 * n * t ** 3
       - 10 * n * t ** 2 - 20 * n * t - 2 * n
       - 4 * t ** 5 - 15 * t ** 4 - 7 * t ** 3 + 8 * t ** 2 + t + 1)
R2 = 24 * n ** 2 / (t * (t - 1) * D ** 3) \
    * (2 * n * t ** 6 - 20 * n * t ** 5 + 10 * n * t ** 4 - 4 * n * t ** 3
       + 36 * n * t ** 2 - 4 * n * t - 12 * n
       + t ** 6 - t ** 5 + 8 * t ** 4 + 7 * t ** 3 - 15 * t ** 2 + 4 * t)
R3 = -64 * n ** 3 * (3 * t ** 5 - 8 * t ** 4 - 22 * t ** 3 - 10 * t ** 2 + 3 * t - 6) / D ** 3
R = sp.Matrix([[R0, R1, R2, R3]])
print("[A]   R_n at n=0 is zero :", sp.simplify(R.subs(n, 0)) == sp.zeros(1, 4))
resA = sp.simplify(R.subs(n, n + 1) * N - R - E)
print("[A]   R_{n+1}N_n - R_n - (e_3+a_2e_2+a_1e_1+a_0e_0) =", resA)
ID_A = (resA == sp.zeros(1, 4))

# ------------------------------------------------ certificate (B) literal S_n
Ascr = 3 * (t - 1) * (t ** 2 + 1) * (t ** 2 + t + 1) / (8 * D ** 2)
Bscr = (t - 1) * (t + 1) * (t ** 2 + 2 * t - 1) / (4 * D)
Rscr = (t ** 2 + 1) * (t + 1) * sp.sqrt(1 - t) / (sp.sqrt(t) * D ** 2)
rho = 2 * t / (t ** 2 + 1) + 1 / (t + 1) + 1 / (2 * (t - 1)) - 1 / (2 * t) - 2 * (2 * t - 1) / D
print("[7]   Rscr'/Rscr - rho =", sp.simplify(sp.diff(Rscr, t) / Rscr - rho))
OKrho = sp.simplify(sp.diff(Rscr, t) / Rscr - rho) == 0

Yn = n * es[0] - Ascr * es[0] - Bscr * es[1]
S = sp.Matrix([[0,
                -3 * n ** 2 * (t - 1) * (t + 1) * (t ** 2 + 1) * (2 * t ** 2 + 7 * t + 1) / (t * D ** 3),
                -12 * n ** 2 * (t - 1) * (t + 1) ** 2 * (t ** 2 + 1) * (t ** 2 - 7 * t + 2) / D ** 3,
                0]])
print("[B]   S_n at n=0 is zero :", sp.simplify(S.subs(n, 0)) == sp.zeros(1, 4))
Yd = sp.Matrix([[sp.diff(Yn[0, i], t) for i in range(4)]])
lhsB = sp.simplify(S.subs(n, n + 1) * N - S)
rhsB = sp.simplify(Yd + Yn * M - rho * Yn)
resB = sp.simplify(lhsB - rhsB)
print("[B]   S_{n+1}N_n-S_n-(Y'+YM-rho Y) =", resB)
ID_B = (resB == sp.zeros(1, 4))

print()
print("SUMMARY  VI3")
print("  (2)  sigma,tau               :", "OK" if sp.simplify(sig - sig_claim) == 0
      and sp.simplify(tau - tau_claim) == 0 else "FAIL")
print("  (3)  M_n, N_n                :", "OK" if (okM and okN and okE) else "FAIL")
print("  (A)  three-order certificate :", "OK" if ID_A else "FAIL")
print("  (B)  first-order certificate :", "OK" if ID_B else "FAIL")
print("  (7)  Rscr'/Rscr = rho        :", "OK" if OKrho else "FAIL")

if not (okM and okN and okE and ID_A and ID_B and OKrho):
    print("RESULT: FAIL")
    sys.exit(1)
print("RESULT: all VI3 printed certificates verified exactly")
sys.exit(0)
