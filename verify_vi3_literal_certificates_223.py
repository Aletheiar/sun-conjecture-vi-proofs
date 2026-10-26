"""Verify the human-readable VI3 certificate appendix without reading JSON.

The formulas below are literal transcriptions of the appendix, not extracted
from the construction outputs. All checks are rational-function identities.
"""
import sympy as s
t,n,v=s.symbols("t n v")
D=t*t-t-1
b=2*(t*t+t-1); d=16*t*(t*t-1)
w=1/(64*(t-1)*(t+1)**2)
sigma=(2*t-1)/D; tau=-4*(t*t+1)/D
den=1-2*b*v+d*v*v
assert s.factor((s.diff(b,t)*v-s.diff(d,t)*v*v/2)/den
                -sigma*v*(b-d*v)/den-tau*v*(v*(b-d*v)/den+1))==0
M=s.zeros(4); N=s.zeros(4)
for j in range(4):
    M[j,j]=(3-j)*n*sigma+j*((n-1)*sigma+(2*n-1)*b*tau/d)+n*s.diff(w,t)/w
    if j<3: M[j,j+1]=(3-j)*n*tau
    if j: M[j,j-1]=-j*n*tau/d
    for k in range(4-j):
        N[j,k]=w*s.binomial(3-j,k)*((2*n+1)*b/(n+1))**(3-j-k)*(-n*d/(n+1))**k
def simp(m): return m.applyfunc(s.cancel)
def zero(m): assert all(s.factor(x)==0 for x in m)
e0=s.Matrix([[1,0,0,0]])
rows=[e0]
for j in range(3): rows.append(simp(rows[-1].diff(t)+rows[-1]*M))
a2=3*(2*t-1)*(3*t*t-3*t-1)/(2*t*(t-1)*D)
a1=(75*t**5-3*t**4-128*t**3-t*t+59*t-6)/(4*t*(t-1)**2*(t+1)**2*D)
a0=(60*t**6+21*t**5-130*t**4-35*t**3+65*t*t+50*t-23)/(8*t*(t-1)**3*(t+1)**3*D)
R=s.Matrix([[
 n**3*(6*t**5+3*t**4+10*t**3-22*t*t+8*t+3)/(t*t*(t-1)**2*(t+1)*D**3),
 6*n*n*(12*n*t**6-4*n*t**5-36*n*t**4-4*n*t**3-10*n*t*t-20*n*t-2*n
         -4*t**5-15*t**4-7*t**3+8*t*t+t+1)/(t*t*(t-1)*(t+1)*D**3),
 24*n*n*(2*n*t**6-20*n*t**5+10*n*t**4-4*n*t**3+36*n*t*t-4*n*t-12*n
         +t**6-t**5+8*t**4+7*t**3-15*t*t+4*t)/(t*(t-1)*D**3),
 -64*n**3*(3*t**5-8*t**4-22*t**3-10*t*t+3*t-6)/D**3
]])
zero(R.subs(n,n+1)*N-R-(rows[3]+a2*rows[2]+a1*rows[1]+a0*e0))
zero(R.subs(n,0))
print("E1: literal third-order certificate and zero boundary verified.",flush=True)
AA=3*(t-1)*(t*t+1)*(t*t+t+1)/(8*D**2)
BB=(t-1)*(t+1)*(t*t+2*t-1)/(4*D)
rho=2*t/(t*t+1)+1/(t+1)+1/(2*(t-1))-1/(2*t)-2*(2*t-1)/D
Y=simp(n*e0-AA*e0-BB*rows[1])
S=s.Matrix([[
 0,
 -3*n*n*(t-1)*(t+1)*(t*t+1)*(2*t*t+7*t+1)/(t*D**3),
 -12*n*n*(t-1)*(t+1)**2*(t*t+1)*(t*t-7*t+2)/D**3,
 0
]])
zero(S.subs(n,n+1)*N-S-(Y.diff(t)+Y*M-rho*Y))
zero(S.subs(n,0))
print("E2: literal normal-residual certificate and zero boundary verified.",flush=True)
print("No JSON or guessed-series data were read. Analytic boundary proof remains in the manuscript.",flush=True)
