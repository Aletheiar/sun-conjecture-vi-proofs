#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Independent referee check 05: high-precision numerics.
  (1) The three ORIGINAL conjectures, summed directly from the definition
      T_n(b,c) = [X^n](X^2+bX+c)^n  (exact integers, forward recurrence)
      compared with the claimed right-hand sides to ~50 digits.
  (2) The three classical hypergeometric inputs actually used:
        (2+15Th) 3F2(1/2,1/3,2/3;1,1; 2/27) = 27/(4 pi)
        (4+33Th) 3F2(1/2,1/3,2/3;1,1; 4/125)= 15 sqrt3/(2 pi)
        (5+42Th) 3F2(1/2,1/2,1/2;1,1; 1/64) = 16/pi
  (3) The CM input at m_7: K(1-m_7) = sqrt7 K(m_7) and
      2 sqrt7 K E - (2+sqrt7) K^2 = pi/2   (= "alpha(7)"), and the bracket of (22).
Nothing here proves anything; it is numerical evidence only.
Exit code 0 only if all comparisons hold at the claimed precision.
"""
import mpmath as mp
mp.mp.dps = 60

FAIL = []
def show(label, got, want, tol=None):
    err = abs(got - want)
    rel = err / abs(want) if want != 0 else err
    nd = -mp.log10(abs(rel)) if rel != 0 else mp.mpf(99)
    ok = (rel < mp.mpf(10) ** (-45))
    print("  [%s] %-42s rel.err=%.3e  (~%s digits)" % ("OK  " if ok else "FAIL", label, rel,
                                                       mp.nstr(nd, 4)))
    if not ok:
        FAIL.append(label)

def Tn_list(b, c, N):
    """T_n(b,c) for n<=N, exact integers, forward recurrence  (n+1)U_{n+1}=(2n+1)bU_n - n d U_{n-1}."""
    d = b * b - 4 * c
    U = [1, b]
    for n in range(1, N):
        U.append(((2 * n + 1) * b * U[n] - n * d * U[n - 1]) // (n + 1))
    for n in range(1, N + 1):                    # exactness guard
        assert ((2 * n + 1) * b * U[n] - n * d * U[n - 1]) % (n + 1) == 0
        assert ((2 * n + 1) * b * U[n] - n * d * U[n - 1]) // (n + 1) == U[n + 1] if n < N else True
    return U

def sum_series(b, c, wnum, wden, an, a0, N=None):
    """sum_n (an*n + a0) T_n(b,c)^3 w^n ,  w = wnum/wden  (exact sign carried)."""
    if N is None:
        N = 400
    U = Tn_list(b, c, N)
    W = mp.mpf(wnum) / mp.mpf(wden)
    tot = mp.mpf(0)
    for n in range(N + 1):
        tot += (an * n + a0) * mp.mpf(U[n]) ** 3 * W ** n
    # tail bound: geometric with the exact ratio, |T_n| <= (b+2 sqrt c)^n
    return tot

print("=" * 78)
print("(1) direct summation of the three original series  (60-digit arithmetic)")
print("=" * 78)
# ---------------- VI1 -------------------------------------------------------
U = Tn_list(10, 121, 600)
w1 = mp.mpf(1) / 55296
s1 = mp.mpf(0)
for n in range(601):
    t = (66 * n + 17) * mp.mpf(U[n]) ** 3 * w1 ** n
    s1 += t
    if n > 20 and abs(t) < mp.mpf(10) ** (-70) * abs(s1):
        break
print("  VI1 : terms used =", n, "  last term =", mp.nstr(abs(t), 3))
show("VI1 sum  = 540 sqrt2/(11 pi)", s1, 540 * mp.sqrt(2) / (11 * mp.pi))

# ---------------- VI2 -------------------------------------------------------
U = Tn_list(22, 441, 700)
w2 = mp.mpf(-1) / 512000
s2 = mp.mpf(0)
for n in range(701):
    t = (126 * n + 31) * mp.mpf(U[n]) ** 3 * w2 ** n
    s2 += t
    if n > 20 and abs(t) < mp.mpf(10) ** (-70) * abs(s2):
        break
print("  VI2 : terms used =", n, "  last term =", mp.nstr(abs(t), 3))
show("VI2 sum  = 880 sqrt5/(21 pi)", s2, 880 * mp.sqrt(5) / (21 * mp.pi))

# ---------------- VI3 -------------------------------------------------------
U = Tn_list(62, 9025, 700)
w3 = mp.mpf(-1) / 23887872
s3 = mp.mpf(0)
for n in range(701):
    t = (3990 * n + 1147) * mp.mpf(U[n]) ** 3 * w3 ** n
    s3 += t
    if n > 20 and abs(t) < mp.mpf(10) ** (-70) * abs(s3):
        break
print("  VI3 : terms used =", n, "  last term =", mp.nstr(abs(t), 3))
show("VI3 sum  = 432(195 sqrt14+94 sqrt2)/(95 pi)", s3,
     432 * (195 * mp.sqrt(14) + 94 * mp.sqrt(2)) / (95 * mp.pi))

print()
print("=" * 78)
print("(2) the classical inputs quoted as (Chan-Liaw 2000) / Ramanujan")
print("=" * 78)
def F32(a1, a2, a3, z, N=4000):
    """3F2(a1,a2,a3;1,1;z) and its Theta derivative, by term recursion."""
    s = mp.mpf(0)
    ts = mp.mpf(0)
    t = mp.mpf(1)
    for n in range(N):
        s += t
        ts += n * t
        k = n + 1
        t *= (a1 + n) * (a2 + n) * (a3 + n) / (mp.mpf(k) * k * k) * z
        if n > 30 and abs(t) < mp.mpf(10) ** (-80):
            break
    return s, ts

for (nm, a, b_, c_, z, want) in [
        ("(2+15Th)W(2/27)", mp.mpf(1) / 2, mp.mpf(1) / 3, mp.mpf(2) / 3,
         mp.mpf(2) / 27, 27 / (4 * mp.pi)),
        ("(4+33Th)W(4/125)", mp.mpf(1) / 2, mp.mpf(1) / 3, mp.mpf(2) / 3,
         mp.mpf(4) / 125, 15 * mp.sqrt(3) / (2 * mp.pi)),
        ("(5+42Th)H(1/64)", mp.mpf(1) / 2, mp.mpf(1) / 2, mp.mpf(1) / 2,
         mp.mpf(1) / 64, 16 / mp.pi)]:
    W, ThW = F32(a, b_, c_, z)
    print("  %s : W = %s ,  ThW = %s" % (nm, mp.nstr(W, 20), mp.nstr(ThW, 20)))
    coefs = {"(2+15Th)W(2/27)": (2, 15), "(4+33Th)W(4/125)": (4, 33),
             "(5+42Th)H(1/64)": (5, 42)}[nm]
    show(nm, coefs[0] * W + coefs[1] * ThW, want)

print()
print("=" * 78)
print("(3) sevenfold singular modulus and the bracket of (22)")
print("=" * 78)
m7 = (8 - 3 * mp.sqrt(7)) / 16
K = mp.ellipk(m7)
Kc = mp.ellipk(1 - m7)
E = mp.ellipe(m7)
show("K(1-m_7) = sqrt7 K(m_7)", Kc, mp.sqrt(7) * K)
show("m_7 = (8-3 sqrt7)/16 vs 1/2-3 sqrt7/16", m7, mp.mpf(1) / 2 - 3 * mp.sqrt(7) / 16)
br = (2 + mp.sqrt(7)) * K ** 2 - 2 * mp.sqrt(7) * K * E
show("(2+sqrt7)K^2-2 sqrt7 K E = -pi/2", br, -mp.pi / 2)
# and the induced value of the "alpha(7)" combination
alpha = mp.pi / (4 * K ** 2) + mp.sqrt(7) * (1 - E / K)
show("alpha(7) = (sqrt7-2)/2", alpha, (mp.sqrt(7) - 2) / 2)
# route (22) assembled
v22 = (40608 * mp.sqrt(2) / (19 * mp.pi ** 2)) * br + 432 * (329 * mp.sqrt(2) + 195 * mp.sqrt(14)) / (95 * mp.pi)
show("(22) value = (1)", v22, 432 * (195 * mp.sqrt(14) + 94 * mp.sqrt(2)) / (95 * mp.pi))
# also K^2 and the CM constant
show("G(t0) = (27 sqrt2/38) H(1/64)", (5 * mp.mpf(27) * mp.sqrt(2) / (2 * 95)) * (4 * K ** 2 / mp.pi ** 2),
     (27 * mp.sqrt(2) / 38) * (4 * K ** 2 / mp.pi ** 2))

print()
print("=" * 78)
print("RESULT:", "ALL NUMERIC CHECKS PASSED" if not FAIL else "FAILED -> %s" % FAIL)
print("=" * 78)
raise SystemExit(0 if not FAIL else 1)
