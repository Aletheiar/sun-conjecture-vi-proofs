#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Independent referee check, script 03 -- do the two operators really annihilate
the elliptic / hypergeometric expressions?

PART A (VI1/VI2)
  claim (2):  G = P W(z),  z = 108 x^2 (1-4x),  P = sqrt(1-4x)/(1-3x),
              W = _3F_2(1/2,1/3,2/3 ; 1,1 ; z),
              theta_z^3 W = z (theta_z+1/2)(theta_z+1/3)(theta_z+2/3) W.
  Check  L(PW) = 0  identically, L = eq (6), by reducing modulo the W-equation.
  P enters only through p'/p (rational), so the computation is rational.

PART B (VI3)
  claim (10): G = J = (P/pi^2)(3K(m)^2 + K(1-m)^2),  m = (1-sqrt t)^2/(2(1+t)),
              P = sqrt(1-t)/(1+t-t^2), K = complete elliptic integral of the
              first kind in the PARAMETER m.
  Check  L(J) = 0  identically, L = eq (5), reducing modulo the two Legendre
  equations satisfied by K(m) and K(1-m).  u = sqrt t makes coefficients rational.
"""
import sys
import sympy as sp

x = sp.symbols('x')
FAIL = []

def check(tag, expr):
    val = sp.simplify(sp.together(sp.expand(expr)))
    ok = (val == 0)
    print(("  %-46s" % tag) + (":  0" if ok else ":  NOT ZERO -> " + str(sp.factor(val))))
    if not ok:
        FAIL.append(tag)
    return ok

def stirling2(k, j):
    if j < 0 or j > k:
        return 0
    if k == 0 and j == 0:
        return 1
    return j * stirling2(k - 1, j) + stirling2(k - 1, j - 1)

# =========================================================== PART A : VI1/VI2
print("=" * 78)
print("PART A : VI1/VI2,  L(PW) = 0 ?")
print("=" * 78)

c3 = 1 - 3 * x - 96 * x ** 2 + 800 * x ** 3 - 2304 * x ** 4 + 2304 * x ** 5
c2 = 3 * x - 276 * x ** 2 + 2928 * x ** 3 - 10944 * x ** 4 + 13824 * x ** 5
c1 = -x - 276 * x ** 2 + 3964 * x ** 3 - 18144 * x ** 4 + 27072 * x ** 5
c0 = -x - 102 * x ** 2 + 1944 * x ** 3 - 10368 * x ** 4 + 17280 * x ** 5

z = 108 * x ** 2 * (1 - 4 * x)
rho = sp.simplify(z / sp.diff(z, x))
print("z =", z, "   rho = z/z' =", rho)

s1 = sp.Rational(1, 2) + sp.Rational(1, 3) + sp.Rational(2, 3)
s2 = sp.Rational(1, 6) + sp.Rational(1, 3) + sp.Rational(2, 9)
s3 = sp.Rational(1, 9)
print("s1,s2,s3 =", s1, s2, s3)

rp, rpp = sp.diff(rho, x), sp.diff(rho, x, 2)
# (1-z) theta_z^3 W = z (s1 theta_z^2 + s2 theta_z + s3) W
A3 = sp.expand((1 - z) * rho ** 3)
A2 = sp.expand((1 - z) * 3 * rho ** 2 * rp - z * s1 * rho ** 2)
A1 = sp.expand((1 - z) * (rho * rho * rpp + rho * rp ** 2) - z * s1 * rho * rp - z * s2 * rho)
A0 = sp.expand(-z * s3)
g2, g1, g0 = [sp.simplify(-A / A3) for A in (A2, A1, A0)]
print("W''' = g2 W'' + g1 W' + g0 W ;  g2 =", g2)

Gm = sp.Matrix([[0, 1, 0], [0, 0, 1], [g0, g1, g2]])
vs = [sp.Matrix([[1, 0, 0]])]
for j in range(3):
    v = vs[-1]
    vd = sp.Matrix([[sp.diff(v[0, i], x) for i in range(3)]])
    vs.append(sp.simplify(vd + v * Gm))

Pexpr = sp.sqrt(1 - 4 * x) / (1 - 3 * x)
pi1 = sp.simplify(sp.log(Pexpr).diff(x))
print("p'/p =", pi1)
pis = [sp.Integer(1)]
for j in range(3):
    pis.append(sp.simplify(sp.diff(pis[-1], x) + pis[-1] * pi1))

ell = [sp.expand(c3 * stirling2(3, j) * x ** j + c2 * stirling2(2, j) * x ** j
                 + c1 * stirling2(1, j) * x ** j + c0 * stirling2(0, j) * x ** j)
       for j in range(4)]
res = sp.zeros(1, 3)
for j in range(4):
    for i in range(j + 1):
        res += ell[j] * sp.binomial(j, i) * pis[j - i] * vs[i]
okA = True
for i in range(3):
    okA &= check("coefficient of D_x^%d W" % i, res[0, i])
print("PART A :", "L(PW) = 0 verified exactly" if okA else "FAILED")

# =============================================================== PART B : VI3
print()
print("=" * 78)
print("PART B : VI3,  3K(m)^2 + K(1-m)^2  annihilated by L ?")
print("=" * 78)

u = sp.symbols('u', positive=True)
tt = u ** 2                      # tt plays the role of t
D_ = tt ** 2 - tt - 1
m = (1 - u) ** 2 / (2 * (1 + tt))
print("m(t) =", sp.simplify(m), "  (t = u^2)")
mu, muu = sp.diff(m, u), sp.diff(m, u, 2)

# m(1-m) k'' + (1-2m) k' - k/4 = 0 ;  for y = K(m):  y_{uu} = A1 y_u + B y
#                                      for y = K(1-m): y_{uu} = A2 y_u + B y
Bcoef = sp.simplify(mu ** 2 / (4 * m * (1 - m)))
A1coef = sp.simplify(muu / mu + mu * (2 * m - 1) / (m * (1 - m)))
A2coef = sp.simplify(muu / mu + mu * (2 * m - 1) / (m * (1 - m)))   # K(1-m) solves the SAME 2nd-order ODE -> A2 = A1
print("B  =", sp.simplify(Bcoef))
print("A1 =", sp.simplify(A1coef))
print("A2 =", sp.simplify(A2coef))

a2 = 3 * (2 * tt - 1) * (3 * tt ** 2 - 3 * tt - 1) / (2 * tt * (tt - 1) * D_)
a1 = (75 * tt ** 5 - 3 * tt ** 4 - 128 * tt ** 3 - tt ** 2 + 59 * tt - 6) \
    / (4 * tt * (tt - 1) ** 2 * (tt + 1) ** 2 * D_)
a0 = (60 * tt ** 6 + 21 * tt ** 5 - 130 * tt ** 4 - 35 * tt ** 3 + 65 * tt ** 2 + 50 * tt - 23) \
    / (8 * tt * (tt - 1) ** 3 * (tt + 1) ** 3 * D_)

def op_mul_left_alphaD(A, alpha, uu):
    """alpha*D_u  composed (on the left) with A = sum A[k] D_u^k."""
    out = [sp.Integer(0)] * (len(A) + 1)
    for k, ak in enumerate(A):
        if ak == 0:
            continue
        out[k] += alpha * sp.diff(ak, uu)
        out[k + 1] += alpha * ak
    return out

def Dt_power(j, alpha, uu):
    A = [sp.Integer(1)]
    for _ in range(j):
        A = op_mul_left_alphaD(A, alpha, uu)
    return [sp.simplify(a) for a in A]

alpha = 1 / (2 * u)
DtP = [Dt_power(j, alpha, u) for j in range(4)]
for j in range(4):
    print("   D_t^%d = sum_k %s D_u^k" % (j, DtP[j]))
# spot checks against D_t f = f'/(2u)
for testf, want in [(u ** 2, 0), (u ** 4, sp.Integer(0)), (u ** 6, sp.Integer(6))]:
    got = sp.simplify(sum(DtP[3][k] * sp.diff(testf, u, k) for k in range(4)))
    print("   spot check D_t^3(%s) = %s   (expected %s)" % (testf, got, want))
    if sp.simplify(got - want) != 0:
        FAIL.append("D_t^3 spot check %s" % testf)

PexprB = sp.sqrt(1 - tt) / (1 + tt - tt ** 2)
pi1B = sp.simplify(sp.log(PexprB).diff(u))     # D_u P = pi1B * P
print("D_u P / P =", sp.simplify(pi1B))
print("   still contains a square root? ", pi1B.has(sp.sqrt(1 - tt)))

pisB = [sp.Integer(1)]
for j in range(3):
    pisB.append(sp.simplify(sp.diff(pisB[-1], u) + pisB[-1] * pi1B))

def annihilates(Aa, Ab, label):
    Gu = sp.Matrix([[0, 1, 1, 0],
                    [Bcoef, Aa, 0, 1],
                    [Bcoef, 0, Ab, 1],
                    [0, Bcoef, Bcoef, Aa + Ab]])
    vsB = [sp.Matrix([[1, 0, 0, 0]])]
    for j in range(3):
        v = vsB[-1]
        vd = sp.Matrix([[sp.diff(v[0, i], u) for i in range(4)]])
        vsB.append(sp.simplify(vd + v * Gu))
    total = sp.zeros(1, 4)
    aj = [a0, a1, a2, sp.Integer(1)]
    for k in range(4):
        coeff = sum(aj[j] * DtP[j][k] for j in range(k, 4))
        coeff = sp.simplify(coeff)
        for i in range(k + 1):
            total += coeff * sp.binomial(k, i) * pisB[k - i] * vsB[i]
    ok = True
    for i in range(4):
        ok &= check("%s : state component %d" % (label, i), total[0, i])
    return ok

ok1 = annihilates(A1coef, A1coef, "K^2      ")
ok2 = annihilates(A1coef, A2coef, "K K_c    ")
ok3 = annihilates(A2coef, A2coef, "K_c^2    ")
okB = ok1 and ok2 and ok3
print("PART B :", "L(3K^2+K_c^2) = 0 verified exactly" if okB else "FAILED")

print()
if FAIL:
    print("RESULT: FAIL --", FAIL)
    sys.exit(1)
print("RESULT: both annihilation identities verified exactly")
sys.exit(0)
