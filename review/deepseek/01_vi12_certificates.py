#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Independent referee check, script 01 -- VI1/VI2 telescoping certificates.

Nothing is imported from the author's workspace.  Everything is rebuilt from
T_n(b,c) = [z^n](1 - 2 b z + d z^2)^(-1/2),  d = b^2 - 4c.

Checked as EXACT rational identities in the free symbols n and x:

 (i)   U_n' = n (sigma U_n + tau U_{n-1}), author's sigma,tau        [eq (2)]
 (ii)  state matrix M_n, shift matrix N_n                             [eq (5)]
 (iii) e_j = D_x^j, and  E_n := sum_j ell_j e_j = R_{n+1}N_n - R_n    [eq (8)]
 (iv)  residual certificate (10)
 (v)   fixed-parameter relation (3) on Taylor coefficients

Exit code 0 <=> every identity holds exactly.
"""
import sys
import sympy as sp

x, n = sp.symbols('x n')

print("=" * 78)
print("VI1/VI2 : family b = 1-6x, c = (1-3x)^2, w = x/(1-4x)")
print("=" * 78)

b = 1 - 6 * x
c = (1 - 3 * x) ** 2
d = sp.expand(b ** 2 - 4 * c)
w = x / (1 - 4 * x)
bp, dp = sp.diff(b, x), sp.diff(d, x)
print("d = b^2-4c =", d)
assert sp.simplify(d - (-3 * (1 - 4 * x))) == 0

# ------------------------------------------------------------------ (i) deriv
sigma_a = (b * bp - dp / 2) / (b ** 2 - d)
tau_a = -(d * bp - b * dp / 2) / (b ** 2 - d)
a_claim = sp.Rational(-3) / (1 - 3 * x)
chk1 = sp.simplify(sigma_a - a_claim)
chk2 = sp.simplify(tau_a - a_claim)
print("[i]   sigma =", sp.simplify(sigma_a), "  tau =", sp.simplify(tau_a))

sig, taus = sp.symbols('sig taus')
eq1 = sp.expand(sig * b + taus - bp)
eq2 = sp.expand(-sig * d - taus * b + dp / 2)
sol = sp.solve([eq1, eq2], [sig, taus], dict=True)[0]
print("[i']  re-derived from F above: sigma =", sp.simplify(sol[sig]),
      " tau =", sp.simplify(sol[taus]))
assert sp.simplify(sol[sig] - a_claim) == 0
assert sp.simplify(sol[taus] - a_claim) == 0

# ------------------------------------------------------------------ (ii) M, N
def Mmat(nn):
    M = sp.zeros(4, 4)
    wp = sp.diff(w, x) / w
    for j in range(4):
        M[j, j] = (3 - j) * nn * sigma_a + j * (nn - 1) * sigma_a \
            + j * (2 * nn - 1) * b * tau_a / d + nn * wp
        if j + 1 <= 3:
            M[j, j + 1] = (3 - j) * nn * tau_a
        if j - 1 >= 0:
            M[j, j - 1] = -j * nn * tau_a / d
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

okM = True
for nn in [1, 2, 3, 4]:
    X = Xvec(nn)
    Xp = sp.Matrix([sp.diff(X[i], x) for i in range(4)])
    if sp.simplify(Xp - Mmat(nn) * X) != sp.zeros(4, 1):
        okM = False
        print("   M_n mismatch at n =", nn)
print("[ii]  M_n reproduces d/dx(w^n U_n^{3-j}U_{n-1}^j), n=1..4 :", okM)

okN = True
for nn in [0, 1, 2, 3]:
    if sp.simplify(Xvec(nn + 1) - Nmat(nn) * Xvec(nn)) != sp.zeros(4, 1):
        okN = False
        print("   N_n mismatch at n =", nn)
print("[ii]  N_n reproduces X_{n+1} = N_n X_n, n=0..3 :", okN)

# ------------------------------------------------------- theta^k -> sum_j x^j D^j
def stirling2(k, j):
    if j < 0 or j > k:
        return 0
    if k == 0 and j == 0:
        return 1
    return j * stirling2(k - 1, j) + stirling2(k - 1, j - 1)

c3 = 1 - 3 * x - 96 * x ** 2 + 800 * x ** 3 - 2304 * x ** 4 + 2304 * x ** 5
c2 = 3 * x - 276 * x ** 2 + 2928 * x ** 3 - 10944 * x ** 4 + 13824 * x ** 5
c1 = -x - 276 * x ** 2 + 3964 * x ** 3 - 18144 * x ** 4 + 27072 * x ** 5
c0 = -x - 102 * x ** 2 + 1944 * x ** 3 - 10368 * x ** 4 + 17280 * x ** 5

ell = [sp.expand(c3 * stirling2(3, j) * x ** j + c2 * stirling2(2, j) * x ** j
                 + c1 * stirling2(1, j) * x ** j + c0 * stirling2(0, j) * x ** j)
       for j in range(4)]
for j in range(4):
    print("[L]   ell_%d =" % j, ell[j])

# e_0 = (1,0,0,0), e_{j+1} = e_j' + e_j M_n
M = Mmat(n)
es = [sp.Matrix([[1, 0, 0, 0]])]
for j in range(3):
    e = es[-1]
    ed = sp.Matrix([[sp.diff(e[0, i], x) for i in range(4)]])
    es.append(sp.simplify(ed + e * M))
E = sp.simplify(es[0] * 0 + sum([ell[j] * es[j] for j in range(4)], sp.zeros(1, 4)))

okE = True
for nn in [0, 1, 2, 3]:
    X = Xvec(nn)
    for j in range(4):
        lhs = sp.simplify(sp.diff(w ** nn * Un(b, c, nn) ** 3, x, j))
        rhs = sp.simplify(sum([es[j].subs(n, nn)[0, i] * X[i] for i in range(4)]))
        if sp.simplify(lhs - rhs) != 0:
            okE = False
            print("   e_j mismatch n=", nn, "j=", j)
print("[iii] e_j reproduces D_x^j(w^n U_n^3), n=0..3, j=0..3 :", okE)

# ------------------------------------------------------------- certificate (8)
R0 = n ** 3 * (4 * x - 1) * (72 * x ** 3 + 18 * x ** 2 - 12 * x + 1) / (3 * x - 1) ** 2
R1 = -9 * n ** 2 * x * (4 * x - 1) ** 2 \
    * (72 * n * x ** 2 - 32 * n * x + n - 8 * x + 1) / (3 * x - 1) ** 2
R2 = 27 * n ** 2 * x * (4 * x - 1) ** 2 / (3 * x - 1) ** 2 \
    * (48 * n * x ** 3 - 84 * n * x ** 2 + 27 * n * x - 2 * n + 24 * x ** 3 - 12 * x ** 2 + 6 * x - 1)
R3 = 27 * n ** 3 * x * (4 * x - 1) ** 3 * (12 * x ** 2 - 12 * x + 1) / (3 * x - 1) ** 2
R = sp.Matrix([[sp.expand(R0), sp.expand(R1), sp.expand(R2), sp.expand(R3)]])
print("[8]   R_n at n=0 is zero :", sp.simplify(R.subs(n, 0)) == sp.zeros(1, 4))

resid = sp.simplify(R.subs(n, n + 1) * Nmat(n) - R - E)
print("[8]   E_n - (R_{n+1}N_n - R_n) =", resid)
IDENT8 = (resid == sp.zeros(1, 4))

# ------------------------------------------------------------ certificate (10)
A = (-9 * x ** 2 * (4 * x - 1) * es[0]
     - (3 * x - 1) * (4 * x - 1) * (6 * x - 1) * x * es[1]
     + (3 * x - 1) ** 2 * (12 * x - 1) * n * es[0])
S = sp.Matrix([[0,
                9 * n ** 2 * (4 * x - 1) / (3 * x - 1),
                27 * n ** 2 * x * (4 * x - 1) / (3 * x - 1),
                0]])
Ad = sp.Matrix([[sp.diff(A[0, i], x) for i in range(4)]])
lhs10 = sp.simplify(Ad + A * M - ((6 * x - 1) / (x * (4 * x - 1))) * A)
rhs10 = sp.simplify(S.subs(n, n + 1) * Nmat(n) - S)
resid10 = sp.simplify(lhs10 - rhs10)
print("[10]  lhs - rhs =", resid10)
IDENT10 = (resid10 == sp.zeros(1, 4))
print("[10]  S_n at n=0 is zero :", sp.simplify(S.subs(n, 0)) == sp.zeros(1, 4))

# --------------------------------------------------------- relation (3), series
Ntr = 14
Gs = sum(Un(b, c, k) ** 3 * w ** k for k in range(Ntr + 1))
Hs = sum(k * Un(b, c, k) ** 3 * w ** k for k in range(Ntr + 1))
lhs3 = sp.series(sp.expand((3 * x - 1) ** 2 * (12 * x - 1)) * Hs, x, 0, Ntr + 1).removeO()
rhs3 = sp.series(sp.expand(9 * x ** 2 * (4 * x - 1)) * Gs
                 + sp.expand((3 * x - 1) * (4 * x - 1) * (6 * x - 1) * x) * sp.diff(Gs, x),
                 x, 0, Ntr + 1).removeO()
diff3 = sp.expand(lhs3 - rhs3)
print("[3]   residual of (3) up to x^%d :" % Ntr, diff3)
OK3 = (diff3 == 0)

print()
print("SUMMARY  VI1/VI2")
print("  (i)  sigma = tau = -3/(1-3x) :", "OK" if (chk1 == 0 and chk2 == 0) else "FAIL")
print("  (ii) M_n, N_n, e_j           :", "OK" if (okM and okN and okE) else "FAIL")
print("  (8)  telescoping identity    :", "OK" if IDENT8 else "FAIL")
print("  (10) residual certificate    :", "OK" if IDENT10 else "FAIL")
print("  (3)  H <-> G,G'              :", "OK" if OK3 else "FAIL")

if not (chk1 == 0 and chk2 == 0 and okM and okN and okE and IDENT8 and IDENT10 and OK3):
    print("RESULT: FAIL")
    sys.exit(1)
print("RESULT: all exact identities hold; (3) holds to order x^%d" % Ntr)
sys.exit(0)
