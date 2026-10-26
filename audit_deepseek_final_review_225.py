"""Exact follow-up to the independent DeepSeek report (2026-10-26).

No imports from the author's other checkers or the reviewer's scripts.
These checks support the accompanying analytic arguments; finite expansions
alone are explicitly NOT an all-orders proof of the target series.
"""
import sympy as s

n = s.symbols('n', integer=True, nonnegative=True)
x, r = s.symbols('x r')
count = 0


def zero(label, value):
    global count
    reduced = s.factor(value)
    if reduced != 0:
        raise AssertionError((label, reduced))
    count += 1
    print('[E]', label, '= 0')


ratio_gap = (n+1)/(n+2) - (2*n+1)**2/(4*(n+1)**2)
bound_gap = 1/(n+2) - (2*n+1)**2/(4*(n+1)**3)
zero('P1: original ratio gap has denominator (n+1)^2',
     ratio_gap - (3*n+2)/(4*(n+2)*(n+1)**2))
zero('P1: reviewer bound gap has denominator (n+1)^3',
     bound_gap - (3*n+2)/(4*(n+2)*(n+1)**3))
zero('The two gaps differ by the positive factor n+1',
     ratio_gap - (n+1)*bound_gap)
zero('Wallis induction base C_0^2 = 1', s.binomial(0, 0)**2-1)

# Direct finite coefficient reconstruction, not imported cached coefficients.
b, c, w = 1-6*x, (1-3*x)**2, x/(1-4*x)
G, H = s.S.Zero, s.S.Zero
for j in range(8):
    U = sum(s.binomial(j, 2*k)*s.binomial(2*k, k)*b**(j-2*k)*c**k
            for k in range(j//2+1))
    term = s.series(U**3*w**j, x, 0, 8).removeO()
    G += term
    H += j*term
Y = (3*x-1)**2*(12*x-1)*H - 9*x**2*(4*x-1)*G \
    - (3*x-1)*(4*x-1)*(6*x-1)*x*s.diff(G, x)
print('[finite diagnostic] G =', s.series(G, x, 0, 4))
print('[finite diagnostic] H =', s.series(H, x, 0, 4))
zero('Report -54*x^2 is not the actual Y coefficient', s.expand(Y).coeff(x, 2))
zero('Direct Y coefficients through x^7 vanish (finite diagnostic only)',
     s.series(Y, x, 0, 8).removeO())
Y_truncated = (3*x-1)**2*(12*x-1)*x - 9*x**2*(4*x-1)*(1+x) \
    - (3*x-1)*(4*x-1)*(6*x-1)*x
print('[finite diagnostic] Y using only G=1+x,H=x =',
      s.series(Y_truncated, x, 0, 3))
zero('Exact homogeneous solution log derivative',
     1/x - 2/(1-4*x) - (6*x-1)/(x*(4*x-1)))
zero('VI12 positive-side majorant derivative',
     s.diff(27*x*(1-4*x)**2,x)-27*(1-4*x)*(1-12*x))
zero('VI12 negative-side majorant derivative',
     s.diff(-27*x*(1-4*x)**2,x)+27*(1-4*x)*(1-12*x))
zero('VI12 endpoint majorant at x=1/36',
     (27*x*(1-4*x)**2).subs(x,s.Rational(1,36))-s.Rational(16,27))
zero('VI12 endpoint majorant at x=-1/60',
     (-27*x*(1-4*x)**2).subs(x,-s.Rational(1,60))-s.Rational(64,125))

# At m=1/2, Legendre equation implies K''=K, and the Legendre
# relation implies K*K'=pi/2. A,B denote these two positive values.
A, B = s.symbols('K0 K1', positive=True)
t = r**2
delta = -r/(1+r**2)
K = A+B*delta+A*delta**2/2
Kc = A-B*delta+A*delta**2/2
P = s.series(s.sqrt(1-t)/(1+t-t**2), r, 0, 3).removeO()
raw = [P*K**2, P*K*Kc, P*Kc**2]
for j, item in enumerate(raw):
    zero('Raw elliptic product %s has the same nonzero constant K0^2' % j,
         item.subs(r, 0)-A**2)
basis = [P*(K+Kc)**2, P*(Kc**2-K**2), P*(Kc-K)**2]
leads = [4*A**2, 4*A*B, 4*B**2]
for j, item in enumerate(basis):
    expanded = s.series(item, r, 0, 3).removeO().expand()
    for k in range(j):
        zero('Adapted basis %s lower coefficient %s vanishes' % (j, k),
             expanded.coeff(r, k))
    zero('Adapted basis %s leading coefficient' % j, expanded.coeff(r, j)-leads[j])
matrix = s.Matrix([[s.series(item, r, 0, 3).removeO().expand().coeff(r,k)
                    for item in basis] for k in range(3)])
zero('Adapted basis determinant = 64*K0^3*K1^3, nonzero',
     matrix.det()-64*A**3*B**3)
J = s.series(P*(3*K**2+Kc**2)/s.pi**2, r, 0, 3).removeO().expand()
zero('Actual candidate constant coefficient h', J.coeff(r,0)-4*A**2/s.pi**2)
zero('Candidate sqrt(t) coefficient -2/pi using K0*K1=pi/2',
     J.coeff(r,1).subs(B,s.pi/(2*A))+2/s.pi)
zero('Candidate t coefficient j-h/2', J.coeff(r,2)-(4*B**2-2*A**2)/s.pi**2)

# All-n rational parts of the endpoint majorant, no sampling in n or eta.
zero('Uniform second derivative rational bound for n>=2',
     2-3*n*(n-1)/((2*n-1)*(2*n-3))
     -(5*n-3)*(n-2)/((2*n-1)*(2*n-3)))
eta = s.symbols('eta')
zero('Derivative geometric-series majorant',
     s.diff(1/(1-eta), eta)-1/(1-eta)**2)
zero('Second derivative geometric-series majorant',
     s.diff(1/(1-eta), eta,2)-2/(1-eta)**3)
M1 = 1+(1-s.Rational(1,4))**(-2)
M2 = 2+2*(1-s.Rational(1,4))**(-3)
zero('Uniform Taylor remainder constant = 2342/27',
     6*(M1**2+M2)-s.Rational(2342,27))
q = s.symbols('t', positive=True)
geometric_s = (1-q)**2*(1+q)
kabs = q**3/((1-q)*(1+q)**2)
zero('VI3 s decreases from 1 on (0,1)',
     s.diff(geometric_s,q)-(q-1)*(3*q+1))
zero('VI3 second geometric ratio is increasing on (0,1)',
     s.diff(kabs,q)/kabs-(3-q)/(q*(1-q)*(1+q)))
zero('VI3 second ratio at t=7/9 is 343/512',
     kabs.subs(q,s.Rational(7,9))-s.Rational(343,512))
AA = 3*(q-1)*(q**2+1)*(q**2+q+1)/(8*(q**2-q-1)**2)
BB = (q-1)*(q+1)*(q**2+2*q-1)/(4*(q**2-q-1))
rho_scaled = (q**2+1)*(q+1)*s.sqrt(1-q)/(q**2-q-1)**2
zero('Endpoint A(0)=-3/8', AA.subs(q,0)+s.Rational(3,8))
zero('Endpoint B(0)=-1/4', BB.subs(q,0)+s.Rational(1,4))
zero('Endpoint sqrt(t)*R(t) tends to 1', rho_scaled.subs(q,0)-1)
zero('Nonzero derivative constant from independently determined limits',
     1/s.pi-BB.subs(q,0)*(-1/s.pi)-s.Rational(3,4)/s.pi)

# The incorrect unsigned b bound in the report fails for n=2 near t=0.
t0 = s.Rational(1,10)
b0, a0 = 2*(t0**2+t0-1), 1+t0-t0**2
U2 = b0**2+2*a0**2
assert U2 > (b0+2*a0)**2
count += 1
print('[E] report |U_n| <= (b+2*sqrt(c))^n is false for VI3: n=2,t=1/10')
zero('Correct positive bound at t=1/10 is |b|+2sqrt(c)=4(1-t^2)',
     abs(b0)+2*a0-4*(1-t0**2))
print('EXACT CHECKS:', count, '; FAILURES: 0')
