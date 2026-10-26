"""Independent exact certificate check, gauge and pullback computations.

All identities here use rational functions. CM special values are not asserted.
"""
import json
from pathlib import Path
import sympy as s
from sympy.functions.combinatorial.numbers import stirling
x,n,z,th=s.symbols("x n z th")
root=Path(__file__).parent
def simp(m): return m.applyfunc(s.cancel)
def zero(m): assert all(s.cancel(v)==0 for v in m)
b=1-6*x; d=-3*(1-4*x); w=x/(1-4*x); a=-3/(1-3*x)
M=s.zeros(4); N=s.zeros(4)
for j in range(4):
    M[j,j]=(3-j)*n*a+j*a*(n-1+(2*n-1)*b/d)+n/(x*(1-4*x))
    if j<3: M[j,j+1]=(3-j)*n*a
    if j: M[j,j-1]=-j*n*a/d
    for k in range(4-j):
        N[j,k]=w*s.binomial(3-j,k)*((2*n+1)*b/(n+1))**(3-j-k)*(-n*d/(n+1))**k
M=simp(M); N=simp(N)
# The derivative state is an all-n consequence of the generating function.
tt=s.symbols("tt")
den=1-2*b*tt+d*tt**2
assert s.factor((s.diff(b,x)*tt-s.diff(d,x)*tt**2/2)/den
                -a*(tt*(1+tt)*(b-d*tt)/den+tt))==0
print("E0: all-coefficient generating-function derivative identity: exact.",flush=True)
v=s.Matrix([[1,0,0,0]])
rows=[v]
for j in range(3): rows.append(simp(rows[-1].diff(x)+rows[-1]*M))
Ldata=json.loads((root/"vi12_rank19_operator_206.json").read_text())
L=s.sympify(Ldata["euler"],locals={"x":x,"th":th})
Ec=simp(sum((x**j*sum(L.coeff(th,k)*stirling(k,j,kind=2)
                          for k in range(j,4))*rows[j] for j in range(4)),s.zeros(1,4)))
R=s.Matrix([[s.sympify(v,locals={"x":x,"n":n}) for v in
             json.loads((root/"vi12_telescoping_208.json").read_text())["certificate"]]])
assert json.loads((root/"vi12_telescoping_208.json").read_text())["relation"]=="ode"
zero(R.subs(n,n+1)*N-R-Ec); zero(R.subs(n,0))
print("E1: L G telescoping identity and n=0 boundary: exact.",flush=True)
P=-9*x*x*(4*x-1); Q=-(3*x-1)*(4*x-1)*(6*x-1); S=(3*x-1)**2*(12*x-1)
E=simp(P*v+Q*x*rows[1]+S*n*v)
normal=json.loads((root/"vi12_normal_residual_209.json").read_text())
R2=s.Matrix([[s.sympify(v,locals={"x":x,"n":n}) for v in normal["certificate"]]])
op=[s.sympify(v,locals={"x":x}) for v in normal["operator"]]
zero(R2.subs(n,n+1)*N-R2-op[0]*E-op[1]*(E.diff(x)+E*M))
zero(R2.subs(n,0))
assert s.factor(op[0]+(6*x-1)/(x*(4*x-1)))==0 and op[1]==1
# Residual Y = C*x*sqrt(1-4x); its x coefficient vanishes from U0=U1(0)=1.
assert s.expand(P).coeff(x,1)+Q.subs(x,0)+S.subs(x,0)==0
print("E2: normal residual first-order ODE and zero constant: exact.",flush=True)
# Derive state identities directly from actual polynomial coefficients for several n.
U=[s.Integer(1),b]
for j in range(1,7):
    U.append(s.cancel(((2*j+1)*b*U[j]-j*d*U[j-1])/(j+1)))
for j in range(8):
    prev=U[j-1] if j else s.Integer(0)
    assert s.factor(s.diff(U[j],x)-j*a*(U[j]+prev))==0
    X=s.Matrix([w**j*U[j]**(3-k)*prev**k for k in range(4)])
    zero(X.diff(x)-M.subs(n,j)*X)
    if j<7:
        NX=s.Matrix([w**(j+1)*U[j+1]**(3-k)*U[j]**k for k in range(4)])
        zero(N.subs(n,j)*X-NX)
print("E3: independent polynomial state checks, n=0..7: exact.",flush=True)
p=s.sympify(Ldata["second_p"],locals={"x":x})
q=s.sympify(Ldata["second_q"],locals={"x":x})
l=-s.Rational(1,2)*(-3)/(1-3*x)+s.Rational(1,4)*(-4)/(1-4*x)
pf=s.factor(p+2*l); qf=s.factor(q+p*l+s.diff(l,x)+l*l)
assert s.factor(pf-(1+10*x-120*x*x)/(x*(1-4*x)*(1+12*x)))==0
assert s.factor(qf+24*x/(x*(1-4*x)*(1+12*x)))==0
xx=z/(z-3)**2
dx=s.diff(xx,z); ddx=s.diff(dx,z)
pp=s.factor(pf.subs(x,xx)*dx-ddx/dx)
qq=s.factor(qf.subs(x,xx)*dx*dx)
gz=1-z/3
pp2=s.factor(pp+2*s.diff(gz,z)/gz)
qq2=s.factor(qq+pp*s.diff(gz,z)/gz)
assert s.factor(pp2-(3*z*z-20*z+9)/(z*(z-9)*(z-1)))==0
assert s.factor(qq2-(z-3)/(z*(z-9)*(z-1)))==0
print("E4: two gauge transformations and rational pullback: exact.",flush=True)
AA=9*x*x*(4*x-1)/S
BB=(3*x-1)*(4*x-1)*(6*x-1)/S
out={}
for label,x0,wa,wb in [("VI1",s.Rational(1,36),17,66),("VI2",-s.Rational(1,60),31,126)]:
    av=s.factor((wa+wb*AA).subs(x,x0))
    bv=s.factor((wb*BB).subs(x,x0))
    out[label]={"x":str(x0),"G_coefficient":str(av),"thetaG_coefficient":str(bv)}
    print(label,"a G + b theta G:",av,bv,flush=True)
(root/"vi12_weight_translation_210.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
print("No CM value or VI1--VI3 conjecture is claimed proved.",flush=True)
