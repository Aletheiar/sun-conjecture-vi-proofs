#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Independent referee check 06 (v2): functional identities and asymptotics.
  (A) VI1/VI2, x in the proved interval [-1/60,1/36]:
        G(x) = P W(z),   Theta W analytically,   H = A G + B G'
  (B) VI3:
        G(t) = P(3K^2+Kc^2)/pi^2   for several t in (0,7/9]
        G'(t) from J  =  G'(t) by finite differences of the series
        H = A G + B G' + (3/(4 pi)) R   for several t   [eq (16)]
        (18) (which uses Kc=sqrt7 K and the classical 5+42Theta) only at t0=7/9
  (C) t -> 0+ :  H -> 1/(pi sqrt t),  G' -> -1/(pi sqrt t),
        (G - h + (2/pi)sqrt t)/t -> j - h/2 ;  with enough terms (N ~ 45/t)
  (D) the explicit bound (21) at sample points
  (E) majorants |f_n'| <= (25/9)C_n, |f_n''| <= (182/27)C_n
Exit code 0 only if all pass.
"""
import mpmath as mp
mp.mp.dps = 40

FAIL = []
def show(label, got, want, tol_exp=-25):
    rel = abs(got) if want == 0 else abs(got - want) / abs(want)
    ok = rel < mp.mpf(10) ** tol_exp
    print("  [%s] %-56s rel.err=%s" % ("OK  " if ok else "FAIL", label, mp.nstr(rel, 3)))
    if not ok:
        FAIL.append(label)

def W32(z, deriv=0, N=4000):
    """3F2(1/2,1/3,2/3;1,1;z) and Theta^deriv (Theta = z d/dz)."""
    s = mp.mpf(0); t = mp.mpf(1)
    for n in range(N):
        s += t * (mp.mpf(n) ** deriv)
        k = n + 1
        t *= (mp.mpf(1) / 2 + n) * (mp.mpf(1) / 3 + n) * (mp.mpf(2) / 3 + n) / (mp.mpf(k) ** 3) * z
        if abs(t) < mp.mpf(10) ** (-60) * (abs(s) + 1):
            break
    return s

print("=" * 78)
print("(A) VI1/VI2 : G = P W(z),  H = A G + B G'   (x in [-1/60, 1/36])")
print("=" * 78)
def vi12(x0, N):
    x = mp.mpf(x0)
    d = dict(x=x)
    d['b'] = 1 - 6 * x
    d['c'] = (1 - 3 * x) ** 2
    d['w'] = x / (1 - 4 * x)
    d['z'] = 108 * x ** 2 * (1 - 4 * x)
    d['P'] = mp.sqrt(1 - 4 * x) / (1 - 3 * x)
    d['A'] = 9 * x ** 2 * (4 * x - 1) / ((3 * x - 1) ** 2 * (12 * x - 1))
    d['B'] = x * (4 * x - 1) * (6 * x - 1) / ((3 * x - 1) * (12 * x - 1))
    d['dz'] = 216 * x - 1296 * x ** 2
    d['Pp'] = d['P'] * (3 / (1 - 3 * x) - 2 / (1 - 4 * x))
    dd = d['b'] ** 2 - 4 * d['c']
    U = [mp.mpf(1), d['b']]
    for n in range(1, N):
        U.append(((2 * n + 1) * d['b'] * U[n] - n * dd * U[n - 1]) / (n + 1))
    G = mp.mpf(0); H = mp.mpf(0)
    for n in range(N + 1):
        t = mp.mpf(U[n]) ** 3 * d['w'] ** n
        G += t; H += n * t
    return G, H, d

for x0, N in [(mp.mpf(1) / 36, 400), (mp.mpf(-1) / 60, 400),
              (mp.mpf(1) / 100, 400), (mp.mpf(-1) / 40, 400),
              (mp.mpf(1) / 200, 400)]:
    G, H, d = vi12(x0, N)
    W = W32(d['z']); ThW = W32(d['z'], 1)
    show("x=%s : G = P W(z)" % mp.nstr(d['x'], 6), G, d['P'] * W)
    Gp = d['Pp'] * W + d['P'] * ThW * d['dz'] / d['z']
    show("x=%s : H = A G + B G'" % mp.nstr(d['x'], 6), H, d['A'] * G + d['B'] * Gp)
    hh = mp.mpf(10) ** (-7)
    G1 = vi12(d['x'] + hh, N)[0]; G2 = vi12(d['x'] - hh, N)[0]
    show("x=%s : G' from chain rule vs FD of series" % mp.nstr(d['x'], 6), Gp,
         (G1 - G2) / (2 * hh), tol_exp=-9)

print()
print("=" * 78)
print("(B) VI3 : G = J,  H = A G + B G' + (3/(4 pi))R,  and (18) at t0 = 7/9")
print("=" * 78)
def vi3(t0, N):
    t = mp.mpf(t0); u = mp.sqrt(t)
    D = dict(t=t, u=u)
    D['D'] = t ** 2 - t - 1
    D['b'] = 2 * (t ** 2 + t - 1)
    D['c'] = (1 + t - t ** 2) ** 2
    D['w'] = 1 / (64 * (t - 1) * (t + 1) ** 2)
    D['A'] = 3 * (t - 1) * (t ** 2 + 1) * (t ** 2 + t + 1) / (8 * D['D'] ** 2)
    D['B'] = (t - 1) * (t + 1) * (t ** 2 + 2 * t - 1) / (4 * D['D'])
    D['R'] = (t ** 2 + 1) * (t + 1) * mp.sqrt(1 - t) / (u * D['D'] ** 2)
    D['P'] = mp.sqrt(1 - t) / (1 + t - t ** 2)
    D['m'] = (1 - u) ** 2 / (2 * (1 + t))
    D['z'] = (1 - t) ** 2 / (1 + t) ** 2
    D['Pp'] = D['P'] * (-1 / (2 * (1 - t)) - (1 - 2 * t) / (1 + t - t ** 2))
    D['zpz'] = -2 / (1 - t) - 2 / (1 + t)
    D['dm'] = (u - 1) / (2 * u * (1 + t)) - (1 - u) ** 2 / (2 * (1 + t) ** 2)
    dd = D['b'] ** 2 - 4 * D['c']
    U = [mp.mpf(1), D['b']]
    for n in range(1, N):
        U.append(((2 * n + 1) * D['b'] * U[n] - n * dd * U[n - 1]) / (n + 1))
    G = mp.mpf(0); H = mp.mpf(0)
    for n in range(N + 1):
        term = mp.mpf(U[n]) ** 3 * D['w'] ** n
        G += term; H += n * term
    D['G'] = G; D['H'] = H
    return D

def Jval(d, deriv=0):
    m = d['m']; P = d['P']; dm = d['dm']
    K = mp.ellipk(m); Kc = mp.ellipk(1 - m)
    E = mp.ellipe(m); Ec = mp.ellipe(1 - m)
    Kp = (E - (1 - m) * K) / (2 * m * (1 - m))
    Kcp = (Ec - m * Kc) / (2 * (1 - m) * m)
    Q = 3 * K ** 2 + Kc ** 2
    val = P * Q / mp.pi ** 2
    if deriv == 0:
        return val
    Qp = dm * (6 * K * Kp - 2 * Kc * Kcp)
    return d['Pp'] * Q / mp.pi ** 2 + P * Qp / mp.pi ** 2

for t0, N in [(mp.mpf(7) / 9, 500), (mp.mpf(1) / 2, 500), (mp.mpf(1) / 4, 500),
              (mp.mpf(1) / 10, 500), (mp.mpf(3) / 4, 500)]:
    D = vi3(t0, N)
    show("t=%s : G = P(3K^2+Kc^2)/pi^2" % mp.nstr(D['t'], 6), D['G'], Jval(D))
    GpJ = Jval(D, 1)
    hh = mp.mpf(10) ** (-6)
    G1 = vi3(D['t'] + hh, N)['G']; G2 = vi3(D['t'] - hh, N)['G']
    show("t=%s : G'(from J) = FD of series" % mp.nstr(D['t'], 6), GpJ,
         (G1 - G2) / (2 * hh), tol_exp=-8)
    show("t=%s : H = A G + B G' + (3/(4 pi))R" % mp.nstr(D['t'], 6), D['H'],
         D['A'] * D['G'] + D['B'] * GpJ + 3 * D['R'] / (4 * mp.pi))
    rel = abs(D['H'] - (D['A'] * D['G'] + D['B'] * GpJ)) / abs(D['H'])
    print("        [info] dropping the (3/4pi)R term gives rel.err = %s" % mp.nstr(rel, 4))

D = vi3(mp.mpf(7) / 9, 500)
K = mp.ellipk(D['m'])
hh_ = 4 * K ** 2 / mp.pi ** 2
th = D['z'] * mp.diff(lambda zz: 4 * mp.ellipk((1 - mp.sqrt(1 - zz)) / 2) ** 2 / mp.pi ** 2, D['z'])
Gp18 = (mp.mpf(5) / 2 * D['Pp'] * hh_ + mp.mpf(5) / 2 * D['P'] * D['zpz'] * th
        - D['P'] * mp.sqrt(7) * D['dm'] / (2 * mp.pi * D['m'] * (1 - D['m'])))
show("t=7/9 : (18) equals G'(from J)", Gp18, Jval(D, 1))

print()
print("=" * 78)
print("(C) t -> 0+ asymptotics (N ~ 45/t terms)")
print("=" * 78)
gam = mp.gamma(mp.mpf(1) / 4)
hconst = gam ** 4 / (4 * mp.pi ** 3)
jconst = 16 * mp.pi / gam ** 4
cache = {}
def G_of(t0):
    key = mp.nstr(t0, 20)
    if key not in cache:
        cache[key] = vi3(t0, int(45 / t0) + 50)
    return cache[key]

def rate_check(name, func, limit, ts):
    """Verify: func(t) - limit = c sqrt(t) + O(t)  (i.e. O(sqrt t) rate, and the limit)."""
    print("  %s : limit = %s" % (name, mp.nstr(limit, 14)))
    cs = []
    for t0 in ts:
        dev = func(t0) - limit
        c = dev / mp.sqrt(t0)
        cs.append(c)
        print("    t=%-8s f=%s   dev=%s   dev/sqrt(t)=%s"
              % (mp.nstr(t0, 3), mp.nstr(func(t0), 13), mp.nstr(dev, 6), mp.nstr(c, 8)))
    # deviations must shrink like sqrt(t) (ratio ~ sqrt(10) per factor 10) and /sqrt(t) settle
    ok = True
    for a, b in zip(cs[:-1], cs[1:]):
        if abs(b - a) > abs(a) * mp.mpf('0.35'):
            ok = False
    if abs(cs[-1]) > mp.mpf('1.5'):
        ok = False
    if ok:
        print("    [OK]   %s : deviation = c sqrt(t) + O(t), c stabilised -> limit confirmed" % name)
    else:
        print("    [FAIL] %s" % name)
        FAIL.append(name)

TS = [mp.mpf(10) ** (-2), mp.mpf(10) ** (-3), mp.mpf(10) ** (-4)]
rate_check("sqrt(t) H", lambda t0: mp.sqrt(t0) * G_of(t0)['H'], 1 / mp.pi, TS)
rate_check("sqrt(t) G'", lambda t0: mp.sqrt(t0) * Jval(G_of(t0), 1), -1 / mp.pi, TS)
rate_check("(G-h+(2/pi)sqrt t)/t", lambda t0: ((G_of(t0)['G'] - hconst + 2 * mp.sqrt(t0) / mp.pi) / t0),
           jconst - hconst / 2, TS)
print("  also: H = A G + B G' + (3/(4 pi)) R at a small t (independent series check)")
for t0 in [mp.mpf(10) ** (-2), mp.mpf(10) ** (-3)]:
    D = G_of(t0)
    show("t=%s : H = A G + B G' + (3/(4 pi))R" % mp.nstr(t0, 3), D['H'],
         D['A'] * D['G'] + D['B'] * Jval(D, 1) + 3 * D['R'] / (4 * mp.pi), tol_exp=-15)

print()
print("=" * 78)
print("(D) the explicit bound (21) : |G-H(s)-eta D(s)| <= (2342/27) eta^2 H(s)")
print("=" * 78)
def H32(z, N=6000):
    s = mp.mpf(0); t = mp.mpf(1)
    for n in range(N):
        s += t
        k = n + 1
        t *= (mp.mpf(1) / 2 + n) ** 3 / (mp.mpf(k) ** 3) * z
        if abs(t) < mp.mpf(10) ** (-55):
            break
    return s
def B32(z, N=6000):
    s = mp.mpf(0); t = mp.mpf(1)
    for n in range(N):
        s += t
        k = n + 1
        t *= (n - mp.mpf(1) / 2) * (mp.mpf(1) / 2 + n) ** 2 / (mp.mpf(k) ** 3) * z
        if abs(t) < mp.mpf(10) ** (-55):
            break
    return s
for t0 in [mp.mpf(1) / 5, mp.mpf(1) / 10, mp.mpf(1) / 20, mp.mpf(1) / 50, mp.mpf(1) / 200]:
    D = vi3(t0, int(45 / t0) + 100)
    s = (1 - D['t']) ** 2 * (1 + D['t'])
    eta = D['t'] / (1 - D['t'] ** 2)
    Hs = H32(s)
    Ds = -mp.mpf(3) / 2 * (Hs - B32(s))
    lhs = abs(D['G'] - Hs - eta * Ds)
    rhs = mp.mpf(2342) / 27 * eta ** 2 * Hs
    if lhs < rhs:
        print("  [OK]   t=%s : (21) holds, slack %s" % (mp.nstr(D['t'], 4),
                                                        mp.nstr(rhs / lhs if lhs > 0 else mp.inf, 4)))
    else:
        print("  [FAIL] t=%s : (21) violated (%s > %s)" % (mp.nstr(D['t'], 4),
                                                           mp.nstr(lhs, 5), mp.nstr(rhs, 5)))
        FAIL.append("(21) at t=%s" % mp.nstr(D['t'], 4))

print()
print("=" * 78)
print("(E) majorants |f_n'| <= (25/9) C_n ,  |f_n''| <= (182/27) C_n , eta in [0,1/4]")
print("=" * 78)
def Cn(n):
    return mp.binomial(2 * n, n) / mp.mpf(4) ** n
def fn_derivs(n, eta):
    f = mp.mpf(0); d1 = mp.mpf(0); d2 = mp.mpf(0)
    for j in range(n + 1):
        c = Cn(j) * Cn(n - j)
        if j == 0:
            f += c
        elif j == 1:
            f += c * (-eta); d1 += -c
        else:
            f += c * (-eta) ** j
            d1 += -j * c * (-eta) ** (j - 1)
            d2 += j * (j - 1) * c * (-eta) ** (j - 2)
    return f, d1, d2
worst1 = mp.mpf(0); worst2 = mp.mpf(0)
for n in range(1, 41):
    for e in [mp.mpf(0), mp.mpf(1) / 20, mp.mpf(1) / 10, mp.mpf(3) / 20,
              mp.mpf(1) / 5, mp.mpf(1) / 4]:
        f, d1, d2 = fn_derivs(n, e)
        worst1 = max(worst1, abs(d1) / Cn(n)); worst2 = max(worst2, abs(d2) / Cn(n))
print("  sup |f_n'|/C_n  =", mp.nstr(worst1, 8), " <= M1 =", mp.nstr(mp.mpf(25) / 9, 8))
print("  sup |f_n''|/C_n =", mp.nstr(worst2, 8), " <= M2 =", mp.nstr(mp.mpf(182) / 27, 8))
if worst1 <= mp.mpf(25) / 9:
    print("  [OK]   |f_n'| <= M1 C_n")
else:
    print("  [FAIL] |f_n'| <= M1 C_n"); FAIL.append("M1")
if worst2 <= mp.mpf(182) / 27:
    print("  [OK]   |f_n''| <= M2 C_n")
else:
    print("  [FAIL] |f_n''| <= M2 C_n"); FAIL.append("M2")

print()
print("=" * 78)
print("RESULT:", "ALL FUNCTIONAL CHECKS PASSED" if not FAIL else "FAILED -> %s" % FAIL)
print("=" * 78)
raise SystemExit(0 if not FAIL else 1)
