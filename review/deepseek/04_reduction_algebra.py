#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Independent referee check 04: exact algebra of the final reductions.
  (A) VI1/VI2:  from eq (3)  H = A G + B G',  G = P W(z),  z=108x^2(1-4x),
      verify EXACTLY (as identities in x, then specialised)
          17 G + 66 H  = (80 sqrt2/11 ) (2 + 15 Theta) W        at x = 1/36
          31 G + 126 H = (352 sqrt15/189)(4 + 33 Theta) W       at x = -1/60
      where Theta = z d/dz.  Nothing is assumed about the classical sums.
  (B) VI3: at t0 = 7/9 verify exactly
          G(t0) = (10 P/pi^2) K^2                      (uses Kc = sqrt7 K, (10),(17))
          1147 G + 3990 H = -(1269 sqrt2/19)(5 h + 42 th) + (1/pi)(142128 sqrt2/95 + 16848 sqrt14/19)
      with h = Hscr(1/64), th = Theta Hscr(1/64) treated as FREE symbols
      (this is eq (19)); then with the classical (17) 5h+42th = 16/pi the value is
          (40608 sqrt2 + 84240 sqrt14)/(95 pi) = 432(195 sqrt14 + 94 sqrt2)/(95 pi).
      Also the numeric check of (22) is left to script 05.
Exit code 0 only if every check passes.
"""
import sympy as sp

FAIL = []
def chk(cond, label, value=None):
    if cond:
        print("  [OK]   %s" % label)
    else:
        print("  [FAIL] %s   ->  %s" % (label, value))
        FAIL.append(label)

x = sp.symbols('x', positive=True)
print("=" * 78)
print("(A) VI1 / VI2 : reduction of the weighted sums to Theta W")
print("=" * 78)

# data of the family (VI1_VI2_PROOF eq (1))
b = 1 - 6 * x
c = (1 - 3 * x) ** 2
d = sp.expand(b ** 2 - 4 * c)
w = x / (1 - 4 * x)
z = sp.expand(108 * x ** 2 * (1 - 4 * x))
P = sp.sqrt(1 - 4 * x) / (1 - 3 * x)
print("d =", d, "   w =", w, "   z =", z)

# eq (3) solved for H :  H = A G + B G'
A_ = sp.simplify(9 * x ** 2 * (4 * x - 1) / ((3 * x - 1) ** 2 * (12 * x - 1)))
B_ = sp.simplify((3 * x - 1) * (4 * x - 1) * (6 * x - 1) * x / ((3 * x - 1) ** 2 * (12 * x - 1)))
print("A =", sp.factor(A_))
print("B =", sp.factor(B_))

Pp = sp.simplify(sp.diff(P, x))
zp = sp.simplify(sp.diff(z, x))

def reduce_sum(mG, mH, x0, name):
    """mG*G + mH*H written as  alpha(x)*W + beta(x)*Theta W , specialised at x0."""
    aa = sp.simplify(mG * P + mH * (A_ * P + B_ * Pp))
    bb = sp.simplify(mH * B_ * P * zp / z)
    aa0 = sp.simplify(aa.subs(x, x0))
    bb0 = sp.simplify(bb.subs(x, x0))
    print("\n  %s :  coefficient of W  = %s" % (name, sp.nsimplify(sp.radsimp(aa0))))
    print("  %s :  coefficient of ThW= %s" % (name, sp.nsimplify(sp.radsimp(bb0))))
    return aa0, bb0

# ---- VI1 : x0 = 1/36 ------------------------------------------------------
x1 = sp.Rational(1, 36)
chk(sp.simplify(z.subs(x, x1) - sp.Rational(2, 27)) == 0, "VI1: z(1/36) = 2/27",
    sp.simplify(z.subs(x, x1)))
chk(sp.simplify(P.subs(x, x1) - 8 * sp.sqrt(2) / 11) == 0, "VI1: P(1/36) = 8 sqrt2 / 11",
    sp.simplify(P.subs(x, x1)))
chk(sp.simplify(w.subs(x, x1) - sp.Rational(1, 32)) == 0, "VI1: w(1/36) = 1/32",
    sp.simplify(w.subs(x, x1)))
aa, bb = reduce_sum(sp.Integer(17), sp.Integer(66), x1, "17G+66H at 1/36")
target_a = sp.Rational(80, 11) * sp.sqrt(2) * 2
target_b = sp.Rational(80, 11) * sp.sqrt(2) * 15
chk(sp.simplify(aa - target_a) == 0, "VI1: W-coefficient = (80 sqrt2/11)*2", sp.simplify(aa - target_a))
chk(sp.simplify(bb - target_b) == 0, "VI1: ThW-coefficient = (80 sqrt2/11)*15", sp.simplify(bb - target_b))

# ---- VI2 : x0 = -1/60 -----------------------------------------------------
x2 = sp.Rational(-1, 60)
chk(sp.simplify(z.subs(x, x2) - sp.Rational(4, 125)) == 0, "VI2: z(-1/60) = 4/125",
    sp.simplify(z.subs(x, x2)))
chk(sp.simplify(P.subs(x, x2) - 16 * sp.sqrt(15) / 63) == 0, "VI2: P(-1/60) = 16 sqrt15/63",
    sp.simplify(P.subs(x, x2)))
chk(sp.simplify(w.subs(x, x2) + sp.Rational(1, 64)) == 0, "VI2: w(-1/60) = -1/64",
    sp.simplify(w.subs(x, x2)))
aa2, bb2 = reduce_sum(sp.Integer(31), sp.Integer(126), x2, "31G+126H at -1/60")
t_a2 = sp.Rational(352, 189) * sp.sqrt(15) * 4
t_b2 = sp.Rational(352, 189) * sp.sqrt(15) * 33
chk(sp.simplify(aa2 - t_a2) == 0, "VI2: W-coefficient = (352 sqrt15/189)*4", sp.simplify(aa2 - t_a2))
chk(sp.simplify(bb2 - t_b2) == 0, "VI2: ThW-coefficient = (352 sqrt15/189)*33", sp.simplify(bb2 - t_b2))

print()
print("=" * 78)
print("(B) VI3 : value at t0 = 7/9 and eq (19)")
print("=" * 78)
t = sp.symbols('t', positive=True)
b3 = 2 * (t ** 2 + t - 1)
D3 = t ** 2 - t - 1
# eq (7)
A3 = 3 * (t - 1) * (t ** 2 + 1) * (t ** 2 + t + 1) / (8 * D3 ** 2)
B3 = (t - 1) * (t + 1) * (t ** 2 + 2 * t - 1) / (4 * D3)
R3 = (t ** 2 + 1) * (t + 1) * sp.sqrt(1 - t) / (sp.sqrt(t) * D3 ** 2)
P3 = sp.sqrt(1 - t) / (1 + t - t ** 2)
z3 = (1 - t) ** 2 / (1 + t) ** 2
m3 = (1 - sp.sqrt(t)) ** 2 / (2 * (1 + t))

t0 = sp.Rational(7, 9)
print("P(t0) =", sp.radsimp(sp.simplify(P3.subs(t, t0))))
print("A(t0) =", sp.radsimp(sp.simplify(A3.subs(t, t0))))
print("B(t0) =", sp.radsimp(sp.simplify(B3.subs(t, t0))))
print("R(t0) =", sp.radsimp(sp.simplify(R3.subs(t, t0))))
print("z(t0) =", sp.simplify(z3.subs(t, t0)))
print("m(t0) =", sp.nsimplify(sp.radsimp(sp.simplify(m3.subs(t, t0)))))
chk(sp.simplify(m3.subs(t, t0) - (8 - 3 * sp.sqrt(7)) / 16) == 0,
    "m(t0) = (8-3 sqrt7)/16 = m_7", sp.nsimplify(sp.simplify(m3.subs(t, t0))))
chk(sp.simplify(z3.subs(t, t0) - sp.Rational(1, 64)) == 0, "z(t0) = 1/64")
chk(sp.simplify(P3.subs(t, t0) - 27 * sp.sqrt(2) / 95) == 0, "P(t0) = 27 sqrt2/95",
    sp.simplify(P3.subs(t, t0)))

P0 = sp.simplify(P3.subs(t, t0))
A0 = sp.simplify(A3.subs(t, t0))
B0 = sp.simplify(B3.subs(t, t0))
R0 = sp.simplify(R3.subs(t, t0))
P0p = sp.simplify(sp.diff(P3, t).subs(t, t0))
zz0 = sp.simplify((sp.diff(z3, t) / z3).subs(t, t0))
mm0 = sp.simplify(sp.diff(m3, t).subs(t, t0))
mloc = sp.simplify((m3 * (1 - m3)).subs(t, t0))
print("P'(t0) =", sp.radsimp(P0p), "   z'/z =", sp.radsimp(zz0),
      "   m' =", sp.radsimp(mm0), "   m(1-m) =", sp.radsimp(mloc))

# h = Hscr(1/64), th = Theta Hscr(1/64) as free symbols (only (17) used later)
h, th = sp.symbols('h th')
# (10) at t0 with Kc = sqrt7 K :  G = (P/pi^2)(3K^2+7K^2) = 10 P K^2/pi^2, h = 4K^2/pi^2
G0 = sp.Rational(10, 4) * P0 * h        # = (10 P/pi^2) K^2  with h = 4K^2/pi^2
print("  G(t0) = (5/2) P h ;  (5/2)P =", sp.radsimp(sp.Rational(5, 2) * P0))
G0p = (sp.Rational(5, 2) * P0p * h + sp.Rational(5, 2) * P0 * zz0 * th
       - P0 * sp.sqrt(7) * mm0 / (2 * sp.pi * mloc))
# (16):  H = A G + B G' + (3/(4 pi)) R
H0 = sp.simplify(A0 * G0 + B0 * G0p + 3 * R0 / (4 * sp.pi))
S = sp.expand(sp.simplify(1147 * G0 + 3990 * H0))
S = sp.expand(sp.radsimp(S))
print("\n  1147 G + 3990 H  (h,th free) =")
print("   ", sp.collect(sp.expand(S), [h, th, 1 / sp.pi]))

c_h = sp.simplify(sp.expand(S).coeff(h))
c_th = sp.simplify(sp.expand(S).coeff(th))
rest = sp.simplify(sp.expand(S) - c_h * h - c_th * th)
print("\n  coefficient of h  :", sp.radsimp(c_h))
print("  coefficient of th :", sp.radsimp(c_th))
print("  constant part     :", sp.radsimp(rest))
chk(sp.simplify(c_h + sp.Rational(1269, 19) * 5 * sp.sqrt(2)) == 0,
    "(19): coefficient of h equals -1269 sqrt2/19 * 5", sp.radsimp(c_h))
chk(sp.simplify(c_th + sp.Rational(1269, 19) * 42 * sp.sqrt(2)) == 0,
    "(19): coefficient of th equals -1269 sqrt2/19 * 42", sp.radsimp(c_th))
rest_target = (sp.Rational(142128, 95) * sp.sqrt(2) + sp.Rational(16848, 19) * sp.sqrt(14)) / sp.pi
chk(sp.simplify(rest - rest_target) == 0,
    "(19): constant part = (1/pi)(142128 sqrt2/95 + 16848 sqrt14/19)", sp.radsimp(rest - rest_target))

val = sp.simplify(rest - sp.Rational(1269, 19) * sp.sqrt(2) * 16 / sp.pi)   # 5h+42th = 16/pi
target = (sp.Rational(40608, 95) * sp.sqrt(2) + sp.Rational(84240, 95) * sp.sqrt(14)) / sp.pi
chk(sp.simplify(val - target) == 0, "(1) for VI3: value = (40608 sqrt2 + 84240 sqrt14)/(95 pi)",
    sp.radsimp(val - target))
target2 = 432 * (195 * sp.sqrt(14) + 94 * sp.sqrt(2)) / (95 * sp.pi)
chk(sp.simplify(val - target2) == 0, "(1) for VI3: equals 432(195 sqrt14 + 94 sqrt2)/(95 pi)",
    sp.radsimp(val - target2))

print()
print("=" * 78)
print("RESULT:", "ALL EXACT REDUCTION CHECKS PASSED" if not FAIL else "FAILED -> %s" % FAIL)
print("=" * 78)
raise SystemExit(0 if not FAIL else 1)
