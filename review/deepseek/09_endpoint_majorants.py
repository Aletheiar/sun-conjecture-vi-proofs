#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Independent referee check 09: the §8 endpoint majorant chain of VI3, on which the
o(t) -> O(t^{3/2}) upgrade of (15) and the leading-term derivation of H in §5 rest.

Checked here from scratch (no use of the author's code):
 (A) the closed form  f_n(eta) = sum_j (1/2)_j/j! (1/2)_{n-j}/(n-j)! (-eta)^j
     against the integral representation f_n(eta) = (2/pi) int_0^{pi/2}(cos^2-eta sin^2)^n dtheta,
     and against the three-term recurrence that (20) implies
        (n+1) f_{n+1} = (1-eta)(n+1/2) f_n + eta n f_{n-1} ,  f_0 = 1 , f_1 = (1-eta)/2 ;
 (B) the generating function (20)  sum_n f_n(eta) v^n = 1/sqrt((1-v)(1+eta v)) ;
 (C) the three uniform bounds of §8 for 0 <= eta <= eta_0 = 1/4 :
        |f_n| <= 2 C_n ,  |f_n'| <= M_1 C_n ,  |f_n''| <= M_2 C_n ,
        |(f_n^3)''| <= 12 (M_1^2 + M_2) C_n^3 ,
     with M_1 = 1+(1-eta_0)^{-2} = 25/9 , M_2 = 2+2(1-eta_0)^{-3} = 182/27 , C_n = binom(2n,n)/4^n ;
 (D) the explicit remainder bound (21) on (0,1/5].  H(s) is evaluated in the exact Clausen
     form H(s) = 4 K(m)^2/pi^2 , m = (1-sqrt(1-s))/2 (independently verified elsewhere).
Exit code 0 only if every assertion holds.
"""
import mpmath as mp
mp.mp.dps = 40

FAIL = []
def ok(cond, label, info=""):
    print("  [%s] %-64s %s" % ("OK  " if cond else "FAIL", label, info))
    if not cond:
        FAIL.append(label)

eta0 = mp.mpf(1) / 4
M1 = 1 + 1 / (1 - eta0) ** 2
M2 = 2 + 2 / (1 - eta0) ** 3
h = mp.mpf('1e-18')          # step for finite differences

def C(n):
    return mp.binomial(2 * n, n) / mp.mpf(4) ** n

def f_closed(n, eta):
    """f_n(eta) = sum_{j=0}^n (1/2)_j/j! (1/2)_{n-j}/(n-j)! (-eta)^j"""
    n = int(n)
    s = mp.mpf(0)
    for j in range(n + 1):
        s += mp.rf(mp.mpf(1) / 2, j) / mp.factorial(j) * \
             mp.rf(mp.mpf(1) / 2, n - j) / mp.factorial(n - j) * (-eta) ** j
    return s

def f_rec(eta, N):
    """all f_0..f_N by the three-term recurrence implied by (20)"""
    f = [mp.mpf(1)]
    if N >= 1:
        f.append((1 - eta) / 2)
    for n in range(1, N):
        f.append(((1 - eta) * (n + mp.mpf(1) / 2) * f[n] + eta * n * f[n - 1]) / (n + 1))
    return f

def f_integral(n, eta):
    n = int(n)
    g = lambda th: (mp.cos(th) ** 2 - eta * mp.sin(th) ** 2) ** n
    return 2 / mp.pi * mp.quad(g, [0, mp.pi / 2])

print("eta_0 = 1/4 ,  M_1 = 25/9 = %s ,  M_2 = 182/27 = %s" % (mp.nstr(M1, 10), mp.nstr(M2, 10)))
ok(abs(M1 - mp.mpf(25) / 9) < mp.mpf('1e-30'), "M_1 = 25/9", "")
ok(abs(M2 - mp.mpf(182) / 27) < mp.mpf('1e-30'), "M_2 = 182/27", "")

print()
print("=" * 78)
print("(A) closed form / recurrence / integral representation of f_n")
print("=" * 78)
worst = mp.mpf(0)
for n in [0, 1, 2, 3, 5, 8, 13]:
    for eta in [mp.mpf(0), mp.mpf(1) / 8, mp.mpf(1) / 4, -mp.mpf(1) / 4]:
        worst = max(worst, abs(f_closed(n, eta) - f_integral(n, eta)))
ok(worst < mp.mpf('1e-25'), "closed form = integral rep (n<=13, eta in {0,1/8,1/4,-1/4})",
   "worst diff = %s" % mp.nstr(worst, 5))
worst = mp.mpf(0)
for eta in [mp.mpf(0), mp.mpf(1) / 8, mp.mpf(1) / 4, -mp.mpf('0.7')]:
    fr = f_rec(eta, 60)
    for n in range(0, 61):
        worst = max(worst, abs(fr[n] - f_closed(n, eta)))
ok(worst < mp.mpf('1e-30'), "recurrence from (20) = closed form (n<=60)",
   "worst diff = %s" % mp.nstr(worst, 5))
ok(abs(f_closed(0, mp.mpf('0.3')) - 1) < mp.mpf('1e-30'), "f_0 = 1", "")
ok(abs(f_closed(1, mp.mpf('0.3')) - (mp.mpf(1) / 2 - mp.mpf('0.3') / 2)) < mp.mpf('1e-30'),
   "f_1 = 1/2 - eta/2", "")
w1 = max(abs(f_closed(n, mp.mpf(0)) - C(n)) for n in range(0, 40))
w2 = mp.mpf(0)
for n in range(0, 40):
    d = (f_closed(n, h) - f_closed(n, -h)) / (2 * h)
    w2 = max(w2, abs(d + n * C(n) / (2 * n - 1)))
ok(w1 < mp.mpf('1e-25'), "f_n(0) = C_n for n<40", "worst %s" % mp.nstr(w1, 4))
ok(w2 < mp.mpf('1e-15'), "f_n'(0) = -n C_n/(2n-1) for n<40", "worst %s" % mp.nstr(w2, 4))

print()
print("=" * 78)
print("(B) generating function (20)")
print("=" * 78)
wgen = mp.mpf(0)
for eta in [mp.mpf(1) / 4, mp.mpf(1) / 10, -mp.mpf(1) / 3]:
    for v in [mp.mpf('0.3'), mp.mpf('-0.5'), mp.mpf('0.8')]:
        fr = f_rec(eta, 300)
        lhs = sum(fr[n] * v ** n for n in range(0, 301))
        rhs = 1 / mp.sqrt((1 - v) * (1 + eta * v))
        wgen = max(wgen, abs(lhs - rhs))
ok(wgen < mp.mpf('1e-18'), "sum_n f_n(eta) v^n = 1/sqrt((1-v)(1+eta v)) [300 terms]",
   "worst %s" % mp.nstr(wgen, 4))

print()
print("=" * 78)
print("(C) the three uniform bounds of §8 on 0 <= eta <= 1/4")
print("=" * 78)
etas = [eta0 * i / 40 for i in range(41)]
sup0 = sup1 = sup2 = sup3 = mp.mpf(0)
NMAX = 120
for n in range(0, NMAX + 1):
    cn = C(n)
    for eta in etas:
        f0 = f_closed(n, eta)
        sup0 = max(sup0, abs(f0) / cn)
        if 0 < eta < eta0:
            d1 = (f_closed(n, eta + h) - f_closed(n, eta - h)) / (2 * h)
            d2 = (f_closed(n, eta + h) - 2 * f0 + f_closed(n, eta - h)) / h ** 2
            sup1 = max(sup1, abs(d1) / cn)
            sup2 = max(sup2, abs(d2) / cn)
            g2 = abs(6 * f0 * d1 ** 2 + 3 * f0 ** 2 * d2)
            sup3 = max(sup3, g2 / cn ** 3)
print("   observed sup |f_n|/C_n        = %-12s (bound 2)" % mp.nstr(sup0, 8))
print("   observed sup |f_n'|/C_n       = %-12s (bound M_1 = %s)" % (mp.nstr(sup1, 8), mp.nstr(M1, 8)))
print("   observed sup |f_n''|/C_n      = %-12s (bound M_2 = %s)" % (mp.nstr(sup2, 8), mp.nstr(M2, 8)))
print("   observed sup |(f_n^3)''|/C_n^3 = %-11s (bound 12(M_1^2+M_2) = %s)"
      % (mp.nstr(sup3, 8), mp.nstr(12 * (M1 ** 2 + M2), 8)))
ok(sup0 <= 2, "|f_n| <= 2 C_n", "")
ok(sup1 <= M1, "|f_n'| <= M_1 C_n", "")
ok(sup2 <= M2, "|f_n''| <= M_2 C_n", "")
ok(sup3 <= 12 * (M1 ** 2 + M2), "|(f_n^3)''| <= 12(M_1^2+M_2) C_n^3", "")

print()
print("=" * 78)
print("(D) the explicit remainder bound (21) of §8")
print("=" * 78)
def Hscr_clausen(z):
    """H(z) = 4 K(m)^2/pi^2 , m = (1-sqrt(1-z))/2"""
    m = (1 - mp.sqrt(1 - z)) / 2
    return 4 * mp.ellipk(m) ** 2 / mp.pi ** 2
def Bfun(z):
    """3F2(-1/2,1/2,1/2;1,1;z)"""
    s = mp.mpf(0)
    for n in range(0, 6000):
        term = (mp.rf(-mp.mpf(1) / 2, n) / mp.factorial(n)) * \
               (mp.rf(mp.mpf(1) / 2, n) / mp.factorial(n)) ** 2 * z ** n
        s += term
        if n > 50 and abs(term) < mp.mpf('1e-45') * max(abs(s), mp.mpf(1)):
            break
    return s
def Dfun(z):
    return -mp.mpf(3) / 2 * (Hscr_clausen(z) - Bfun(z))

for t in [mp.mpf(1) / 5, mp.mpf(1) / 10, mp.mpf(1) / 20, mp.mpf(1) / 50, mp.mpf(1) / 200]:
    eta = t / (1 - t ** 2)
    s = (1 - t) ** 2 * (1 + t)
    fr = f_rec(eta, 4000)
    K1 = 4000
    G = sum(fr[n] ** 3 * s ** n for n in range(0, K1 + 1))
    tail = abs(fr[K1] ** 3 * s ** K1 / (1 - s)) + mp.mpf('1e-40')
    lhs = abs(G - Hscr_clausen(s) - eta * Dfun(s))
    rhs = mp.mpf(2342) / 27 * eta ** 2 * Hscr_clausen(s)
    print("   t=%-8s eta=%-9s |G-H-eta D|=%-11s bound=%-11s slack=%-9s (tail est %s)"
          % (mp.nstr(t, 4), mp.nstr(eta, 5), mp.nstr(lhs, 5), mp.nstr(rhs, 5),
             mp.nstr(rhs / lhs if lhs > 0 else mp.inf, 5), mp.nstr(tail, 3)))
    ok(lhs <= rhs, "bound (21) at t=%s" % mp.nstr(t, 4))

print()
print("=" * 78)
print("RESULT:", "ENDPOINT MAJORANT CHAIN VERIFIED" if not FAIL else "FAILED -> %s" % FAIL)
print("=" * 78)
raise SystemExit(0 if not FAIL else 1)
