#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Independent referee check 08: the two EXTERNAL classical inputs.

(A) Chan-Liaw, Pacific J. Math. 192 (2000) 219-238, p.219, (1.1) and (1.2):
      (1.1)  27/(4 pi)    = sum_{m>=0} (2 + 15m) (1/2)_m(1/3)_m(2/3)_m / (m!)^3 (2/27)^m
      (1.2)  15 sqrt3/(2 pi) = sum_{m>=0} (4 + 33m) (1/2)_m(1/3)_m(2/3)_m / (m!)^3 (4/125)^m
    (transcribed from the published text; verified here independently to 60 digits).

(B) Borwein-Borwein, Proc. Indian Acad. Sci. (Math. Sci.) 97 (1987) 53-59, p.57, N=7:
      2 k7 k7' = 1/8,   K(k7')/K(k7) = sqrt7,   alpha(7) = (sqrt7 - 2)/2.
    The paper's K has the MODULUS k as variable; the proof uses the PARAMETER m = k^2.
    We verify the transcription, and then verify that (17)'s third relation
      (5 + 42 Theta) H(1/64) = 16/pi ,  H = 3F2(1/2,1/2,1/2;1,1),  Theta = z d/dz
    really FOLLOWS from m7 and alpha(7) by the Legendre-type derivative formula
      dK/dm = (E - (1-m)K) / (2 m (1-m)),
    i.e. that the proof is not silently assuming the very series it wants to prove.
Exit code 0 only if every assertion holds.
"""
import sympy as sp
import mpmath as mp
mp.mp.dps = 60

FAIL = []
def ok(cond, label, info=""):
    print("  [%s] %-64s %s" % ("OK  " if cond else "FAIL", label, info))
    if not cond:
        FAIL.append(label)

# ---------------------------------------------------------------- (A)
print("=" * 78)
print("(A) Chan-Liaw (1.1), (1.2): the two Ramanujan evaluations")
print("=" * 78)

def threeF2_half_thirds(z, N=4000):
    """sum (1/2)_m (1/3)_m (2/3)_m / (m!)^3 z^m  and its Theta = z d/dz."""
    s = mp.mpf(0); t = mp.mpf(0)
    term = mp.mpf(1)
    for m in range(N):
        t += m * term
        s += term
        term *= (mp.mpf(m) + mp.mpf(1) / 2) * (mp.mpf(m) + mp.mpf(1) / 3) * \
                (mp.mpf(m) + mp.mpf(2) / 3) / (mp.mpf(m) + 1) ** 3 * z
        if abs(term) < mp.mpf('1e-80'):
            break
    return s, t

for z0, A, B, rhs, name in [
        (mp.mpf(2) / 27, mp.mpf(2), mp.mpf(15), 27 / (4 * mp.pi), "(1.1) z=2/27, 2+15m"),
        (mp.mpf(4) / 125, mp.mpf(4), mp.mpf(33), 15 * mp.sqrt(3) / (2 * mp.pi), "(1.2) z=4/125, 4+33m")]:
    S, T = threeF2_half_thirds(z0)
    val = A * S + B * T
    rel = abs(val - rhs) / abs(rhs)
    ok(rel < mp.mpf('1e-50'), "%s : A*W + B*Theta W = rhs" % name,
       "rel.err = %s" % mp.nstr(rel, 6))

print()
print("=" * 78)
print("(B) Borwein-Borwein p.57, N=7 : transcription and the (17) derivation")
print("=" * 78)
sqrt7 = sp.sqrt(7)
M7 = (8 - 3 * sqrt7) / 16                      # parameter m7 = k7^2
k7 = (3 - sqrt7) / (4 * sp.sqrt(2))            # modulus
k7c = sp.sqrt(1 - k7 ** 2)                     # complementary modulus
ok(sp.simplify(2 * k7 * k7c - sp.Rational(1, 8)) == 0,
   "2 k7 k7' = 1/8", "k7 = %s" % sp.nsimplify(sp.simplify(k7)))
ok(sp.simplify(k7 ** 2 - M7) == 0, "m7 = k7^2 = (8-3 sqrt7)/16", str(sp.simplify(k7 ** 2)))
ok(sp.simplify(4 * M7 * (1 - M7) - sp.Rational(1, 64)) == 0,
   "z = 4m(1-m) = 1/64  exact (m = m7)", "")

# full-precision mpmath value of the parameter (no truncation)
m = (mp.mpf(8) - 3 * mp.sqrt(7)) / 16
K = mp.ellipk(m); E = mp.ellipe(m)
alpha = mp.pi / (4 * K ** 2) + mp.sqrt(7) * (1 - E / K)
ok(abs(alpha - (mp.sqrt(7) - 2) / 2) < mp.mpf('1e-50'),
   "pi/(4K^2) + sqrt7(1-E/K) = (sqrt7-2)/2  [p.57 value]",
   "diff = %s" % mp.nstr(alpha - (mp.sqrt(7) - 2) / 2, 5))
ok(abs(K / mp.ellipk(1 - m) - mp.mpf(1) / mp.sqrt(7)) < mp.mpf('1e-50'),
   "K(m7)/K(1-m7) = 1/sqrt7  [p.57: m7 = sqrt7]", "")

# Now the algebra: does (5 + 42Theta)H(1/64) = 16/pi FOLLOW from m7 and alpha(7)?
#  H(z) = 4 K(m(z))^2 / pi^2 with z = 4m(1-m)  (Clausen, m the parameter),
#  Theta = z d/dz = [m(1-m)/(1-2m)] d/dm,   dK/dm = (E-(1-m)K)/(2m(1-m)).
Ks, Es = sp.symbols('Ks Es')          # values at m = m7
# (5+42Theta)H = (4/pi^2) * [ 5K^2 + 42 K (E-(1-m)K)/(1-2m) ]
Phi = 5 * Ks ** 2 + 42 * Ks * (Es - (1 - M7) * Ks) / (1 - 2 * M7)
# alpha(7) = (sqrt7-2)/2  <=>  pi = -2(2+sqrt7)K^2 + 4 sqrt7 K E
pi_expr = -2 * (2 + sqrt7) * Ks ** 2 + 4 * sqrt7 * Ks * Es


# We must show  (4 Phi / pi^2) = 16/pi  <=>  4 Phi - 16 pi = 0,  with pi = pi_expr.
resid = sp.simplify(sp.expand(4 * Phi - 16 * pi_expr))
ok(resid == 0,
   "(5+42Theta)H(1/64)=16/pi  follows from m7 & alpha(7) alone",
   "residual = %s" % resid)
# and numerically, end to end:
val = mp.mpf(4) / mp.pi ** 2 * (5 * K ** 2 + 42 * K * (E - (1 - m) * K) / (1 - 2 * m))
ok(abs(val - 16 / mp.pi) < mp.mpf('1e-50'),
   "numeric end-to-end (5+42Theta)H(1/64) = 16/pi", "rel.err = %s"
   % mp.nstr(abs(val - 16 / mp.pi) / (16 / mp.pi), 5))
# the (5+42Theta) series written as a 3F2 sum (what the target formula needs)
def H3F2(z, N=200):
    s = mp.mpf(0); t = mp.mpf(0); term = mp.mpf(1)
    for nn in range(N):
        t += nn * term; s += term
        term *= ((mp.mpf(nn) + mp.mpf(1) / 2) / (mp.mpf(nn) + 1)) ** 3 * z
    return s, t
S, T = H3F2(mp.mpf(1) / 64)
ok(abs(5 * S + 42 * T - 16 / mp.pi) < mp.mpf('1e-50'),
   "(5+42Theta)H(1/64)=16/pi as a 3F2 series sum", "")
ok(abs(S - 4 * K ** 2 / mp.pi ** 2) < mp.mpf('1e-45'),
   "H(1/64) = 4 K(m7)^2/pi^2  (Clausen, m the parameter)", "")

print()
print("=" * 78)
print("RESULT:", "EXTERNAL INPUTS VERIFIED" if not FAIL else "FAILED -> %s" % FAIL)
print("=" * 78)
raise SystemExit(0 if not FAIL else 1)
