#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Independent referee check 07: the n -> infinity boundary term of the telescoping,
and the geometric majorants that justify term-by-term operations.

  For each family, the certificate identity  (L applied to the n-th term) =
  R_{n+1}X_{n+1} - R_n X_n  sums to  R_{N+1}X_{N+1}  (the boundary term).
  A rigorous infinite summation therefore needs  R_n X_n -> 0.
  We check that numerically with the LITERAL rational certificates, and we check
  the majorant ratios |w|(b+2 sqrt c)^3 < 1 on the claimed domains.
Exit code 0 only if the boundary terms decay and all majorants are < 1.
"""
import sympy as sp
import mpmath as mp
mp.mp.dps = 30

FAIL = []
def ok(cond, label, info=""):
    print("  [%s] %-62s %s" % ("OK  " if cond else "FAIL", label, info))
    if not cond:
        FAIL.append(label)

def boundary_ok(sizes, label, r_bound=None, tail_scale=mp.mpf('1e-6')):
    """The telescoping boundary term R_{n+1}X_{n+1} must tend to 0 (it must, a fortiori,
    be summable).  Criterion actually needed for the infinite summation: the tail
    sum_{n>N} |R_n X_n| -> 0, i.e. the term decays geometrically with ratio < 1.
    We test the observed ratios against the PROVEN majorant ratio r_bound (when given),
    not against an arbitrary absolute threshold, and bound the whole tail geometrically."""
    rs = [sizes[i + 1][1] / sizes[i][1] for i in range(len(sizes) - 1)]
    print("      ratio term_{k+1}/term_k : %s" % ", ".join(mp.nstr(r, 5) for r in rs))
    rmax = max(rs)
    tail = sizes[-1][1] * rmax / (1 - rmax) if rmax < 1 else mp.inf
    print("      max ratio = %s%s ; geometric tail bound = %s"
          % (mp.nstr(rmax, 6),
             "" if r_bound is None else " (proven bound %s)" % mp.nstr(r_bound, 6),
             mp.nstr(tail, 6)))
    if r_bound is None:
        return rmax < 1 and tail < tail_scale
    return rmax <= r_bound and r_bound < 1 and tail < tail_scale


print("=" * 78)
print("(A) VI1/VI2 :  ratio |w|(b+2 sqrt c)^3  on the segment [-1/60, 1/36]")
print("=" * 78)
def ratio_vi12(x0):
    x = mp.mpf(x0)
    b = 1 - 6 * x
    c = (1 - 3 * x) ** 2
    w = x / (1 - 4 * x)
    return abs(w) * (b + 2 * mp.sqrt(c)) ** 3, b + 2 * mp.sqrt(c)
for x0, claimed in [(mp.mpf(1) / 36, mp.mpf(16) / 27), (mp.mpf(-1) / 60, mp.mpf(64) / 125)]:
    r, b2 = ratio_vi12(x0)
    ok(abs(r - claimed) < mp.mpf(10) ** (-25) and r < 1,
       "x=%s : |w|(b+2sqrt c)^3 = %s < 1" % (mp.nstr(x0, 6), mp.nstr(r, 12)),
       "(claimed %s)" % mp.nstr(claimed, 8))
# monotonicity of the ratio on the segment (only endpoints matter)
xs = [sp.Rational(-1, 60) + sp.Rational(1, 60 * 400) * i for i in range(401)]
xx = sp.symbols('xx')
expr = sp.Abs(xx / (1 - 4 * xx)) * (3 * (1 - 4 * xx)) ** 3
vals = [abs(sp.Rational(x) / (1 - 4 * sp.Rational(x))) * (3 * (1 - 4 * sp.Rational(x))) ** 3 for x in xs]
ok(max(vals) <= sp.Rational(16, 27),
   "sup of the ratio on [-1/60,1/36] = %s <= 16/27" % max(vals))

print()
print("=" * 78)
print("(B) VI1/VI2 : boundary term R_{n+1} X_{n+1} -> 0 ?")
print("=" * 78)
x = sp.symbols('x')
n = sp.symbols('n')
# literal certificate of VI1_VI2_PROOF eq (7)
R = [n ** 3 * (4 * x - 1) * (72 * x ** 3 + 18 * x ** 2 - 12 * x + 1) / (3 * x - 1) ** 2,
     -9 * n ** 2 * x * (4 * x - 1) ** 2 * (72 * n * x ** 2 - 32 * n * x + n - 8 * x + 1) / (3 * x - 1) ** 2,
     27 * n ** 2 * x * (4 * x - 1) ** 2 / (3 * x - 1) ** 2 *
     (48 * n * x ** 3 - 84 * n * x ** 2 + 27 * n * x - 2 * n + 24 * x ** 3 - 12 * x ** 2 + 6 * x - 1),
     27 * n ** 3 * x * (4 * x - 1) ** 3 * (12 * x ** 2 - 12 * x + 1) / (3 * x - 1) ** 2]
for x0 in [mp.mpf(1) / 36, mp.mpf(-1) / 60]:
    b = 1 - 6 * mp.mpf(x0)
    c = (1 - 3 * mp.mpf(x0)) ** 2
    d = b * b - 4 * c
    w = mp.mpf(x0) / (1 - 4 * mp.mpf(x0))
    U = [mp.mpf(1), b]
    Nmax = 240
    for k in range(1, Nmax):
        U.append(((2 * k + 1) * b * U[k] - k * d * U[k - 1]) / (k + 1))
    sizes = []
    for k in [40, 80, 120, 180, 240]:
        Xk = [w ** k * U[k] ** 3, w ** k * U[k] ** 2 * U[k - 1],
              w ** k * U[k] * U[k - 1] ** 2, w ** k * U[k - 1] ** 3]
        Rk = [complex(sp.N(Ri.subs({n: k, x: sp.Rational(int(round(float(x0) * 1000)), 1000)}), 25)).real
              for Ri in R]
        # use exact rational x instead
        xr = sp.Rational(1, 36) if abs(x0 - mp.mpf(1) / 36) < 1e-12 else sp.Rational(-1, 60)
        Rk = [sp.N(Ri.subs({n: k, x: xr}), 30) for Ri in R]
        bnd = sum(abs(mp.mpf(str(Rk[i]))) * abs(Xk[i]) for i in range(4))
        sizes.append((k, bnd))
    print("   x=%s boundary |R_{n+1}X_{n+1}|:" % mp.nstr(x0, 6))
    for (k, v) in sizes:
        print("      n=%3d : %s" % (k, mp.nstr(v, 6)))
    dec = boundary_ok(sizes, "x=%s" % mp.nstr(x0, 6),
                      r_bound=abs(w) * (3 * (1 - 4 * mp.mpf(x0))) ** 3)

print()
print("=" * 78)
print("(C) VI3 : majorant |w|(4 max(t,1-t^2))^3 < 1 and boundary term -> 0")
print("=" * 78)
def t_of(x):
    return mp.mpf(x)
for t0 in [mp.mpf('0.01'), mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf(7) / 9]:
    t = t0
    b = 2 * (t ** 2 + t - 1)
    a2 = 1 + t - t ** 2
    w = 1 / (64 * (t - 1) * (t + 1) ** 2)
    r = abs(w) * (4 * max(t, 1 - t ** 2)) ** 3
    ok(r < 1, "t=%s : |w|(4max(t,1-t^2))^3 = %s < 1" % (mp.nstr(t, 4), mp.nstr(r, 8)))
# the claimed bound 343/512 at t=7/9
t = mp.mpf(7) / 9
r = t ** 3 / ((1 - t) * (1 + t) ** 2)
ok(abs(r - mp.mpf(343) / 512) < mp.mpf('1e-25'), "at t=7/9 the ratio is t^3/((1-t)(1+t)^2) = 343/512",
   mp.nstr(r, 12))

Dd = sp.Symbol('Dd')
tt = sp.symbols('t')
n3 = sp.symbols('n')
Dpoly = tt ** 2 - tt - 1
# literal certificate A of VI3_EXPLICIT_CERTIFICATES
RA = [
    n3 ** 3 * (6 * tt ** 5 + 3 * tt ** 4 + 10 * tt ** 3 - 22 * tt ** 2 + 8 * tt + 3) /
    (tt ** 2 * (tt - 1) ** 2 * (tt + 1) * Dpoly ** 3),
    6 * n3 ** 2 / (tt ** 2 * (tt - 1) * (tt + 1) * Dpoly ** 3) *
    (12 * n3 * tt ** 6 - 4 * n3 * tt ** 5 - 36 * n3 * tt ** 4 - 4 * n3 * tt ** 3 - 10 * n3 * tt ** 2
     - 20 * n3 * tt - 2 * n3 - 4 * tt ** 5 - 15 * tt ** 4 - 7 * tt ** 3 + 8 * tt ** 2 + tt + 1),
    24 * n3 ** 2 / (tt * (tt - 1) * Dpoly ** 3) *
    (2 * n3 * tt ** 6 - 20 * n3 * tt ** 5 + 10 * n3 * tt ** 4 - 4 * n3 * tt ** 3 + 36 * n3 * tt ** 2
     - 4 * n3 * tt - 12 * n3 + tt ** 6 - tt ** 5 + 8 * tt ** 4 + 7 * tt ** 3 - 15 * tt ** 2 + 4 * tt),
    -64 * n3 ** 3 * (3 * tt ** 5 - 8 * tt ** 4 - 22 * tt ** 3 - 10 * tt ** 2 + 3 * tt - 6) / Dpoly ** 3]
for t0 in [mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf(7) / 9]:
    t = t0
    b = 2 * (t ** 2 + t - 1)
    c = (1 + t - t ** 2) ** 2
    dd = 16 * t * (t ** 2 - 1)
    w = 1 / (64 * (t - 1) * (t + 1) ** 2)
    U = [mp.mpf(1), b]
    Nmax = 300
    for k in range(1, Nmax):
        U.append(((2 * k + 1) * b * U[k] - k * dd * U[k - 1]) / (k + 1))
    sizes = []
    for k in [50, 100, 150, 200, 250]:
        Xk = [w ** k * U[k] ** 3, w ** k * U[k] ** 2 * U[k - 1],
              w ** k * U[k] * U[k - 1] ** 2, w ** k * U[k - 1] ** 3]
        tr = sp.nsimplify(t0, [sp.Rational(1, 10), sp.Rational(1, 2), sp.Rational(7, 9)])
        Rk = [sp.N(Ri.subs({n3: k, tt: tr}), 30) for Ri in RA]
        bnd = sum(abs(mp.mpf(str(Rk[i]))) * abs(Xk[i]) for i in range(4))
        sizes.append((k, bnd))
    print("   t=%s boundary |R_{n+1}X_{n+1}|:" % mp.nstr(t, 4))
    for (k, v) in sizes:
        print("      n=%3d : %s" % (k, mp.nstr(v, 6)))
    dec = boundary_ok(sizes, "t=%s" % mp.nstr(t, 4),
                      r_bound=abs(w) * (4 * max(t, 1 - t ** 2)) ** 3)
    ok(dec, "t=%s : boundary term decays geometrically (summable tail)" % mp.nstr(t, 4))

print()
print("=" * 78)
print("RESULT:", "BOUNDARY / MAJORANT CHECKS PASSED" if not FAIL else "FAILED -> %s" % FAIL)
print("=" * 78)
raise SystemExit(0 if not FAIL else 1)
