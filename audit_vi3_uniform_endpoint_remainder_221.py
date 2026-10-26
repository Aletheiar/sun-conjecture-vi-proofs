"""Independent algebra for an explicit uniform endpoint Taylor bound.

The inequalities for the integrals and convergent geometric sums are written
out in the proof addendum; this script checks constants and symbolic expansions.
"""
from pathlib import Path
import json
import sympy as s
n=s.symbols("n",integer=True,positive=True)
r,eta,z,v=s.symbols("r eta z v")
K0,K1=s.symbols("K0 K1",positive=True)
# Wallis bound C_n <= 1/sqrt(n+1) follows by induction from this positive gap.
wallis_gap=(n+1)/(n+2)-(2*n+1)**2/(4*(n+1)**2)
assert s.factor(wallis_gap-(3*n+2)/(4*(n+2)*(n+1)**2))==0
# Beta ratio for n>=2, bounded by 2. Equality holds at n=2.
beta_ratio=3*n*(n-1)/((2*n-1)*(2*n-3))
assert s.factor(2-beta_ratio-(5*n-3)*(n-2)/((2*n-1)*(2*n-3)))==0
eta0=s.Rational(1,4)
M1=1+(1-eta0)**-2
M2=2+2*(1-eta0)**-3
C=s.factor(6*(M1*M1+M2))
assert (M1,M2,C)==(s.Rational(25,9),s.Rational(182,27),s.Rational(2342,27))
print("E1: Wallis induction, second-derivative beta bound, remainder constant",C,flush=True)
# Independent finite polynomial checks of the exact generating function:
# sum f_n(eta)*v^n = 1/sqrt((1-v)(1+eta*v)).
F=s.series(((1-v)*(1+eta*v))**(-s.Rational(1,2)),v,0,8).removeO()
for j in range(8):
    fj=s.expand(F).coeff(v,j)
    convolution=sum(s.rf(s.Rational(1,2),k)*s.rf(s.Rational(1,2),j-k)
                    *(-eta)**k/(s.factorial(k)*s.factorial(j-k)) for k in range(j+1))
    assert s.expand(fj-convolution)==0
    Cj=s.binomial(2*j,j)/4**j
    assert fj.subs(eta,0)==Cj
    if j: assert s.diff(fj,eta).subs(eta,0)==-s.Rational(j,2*j-1)*Cj
print("E2: independent finite generating-polynomial checks n=0..7 (not an all-n proof).",flush=True)
# Taylor-expand the actual candidate and the comparison hypergeometric branch
# in r=sqrt(t). K''(1/2)=K0 follows directly from the Gauss equation.
t=r*r
P=s.sqrt(1-t)/(1+t-t*t)
delta=-r/(1+r*r)
Km=K0+K1*delta+K0*delta*delta/2
Kc=K0-K1*delta+K0*delta*delta/2
J=s.series(P*(3*Km*Km+Kc*Kc)/s.pi**2,r,0,3).removeO().expand()
ss=(1-t)**2*(1+t)
delta_s=-r*s.sqrt(1+r*r-r**4)/2
Ks=K0+K1*delta_s+K0*delta_s*delta_s/2
Hs=s.series(4*Ks*Ks/s.pi**2,r,0,3).removeO().expand()
h=4*K0*K0/s.pi**2
jc=4*K1*K1/s.pi**2
D1=-3*(h-jc)/4
expected=s.expand(Hs+D1*r*r)
assert s.simplify(J-expected)==0
assert J.coeff(r,0)==h
assert s.simplify(J.coeff(r,1).subs(K1,s.pi/(2*K0))+2/s.pi)==0
assert s.simplify(J.coeff(r,2)-(jc-h/2))==0
print("E3: independently expanded coefficients of 1,sqrt(t),t match exactly.",flush=True)
# Contiguous relation D(s) = -3/2*(H(s)-B(s)), coefficient by coefficient.
an=s.symbols("an")
bn=-an/(2*n-1)
assert s.factor(-3*n*an/(2*n-1)+s.Rational(3,2)*(an-bn))==0
print("E4: all-index rational contiguous coefficient identity: exact.",flush=True)
out={"eta0":"1/4","t_valid":"0<t<=1/5","M1":str(M1),"M2":str(M2),
     "bound_constant":str(C),
     "bound":"abs(G-Hclassic(s)-eta*D(s)) <= (2342/27)*eta^2*Hclassic(s)",
     "consequence":"G=h-2/pi*sqrt(t)+(j-h/2)*t+O(t^(3/2))",
     "scope":"Analytic inequalities in addendum; finite polynomial checks labelled separately"}
Path(__file__).with_name("vi3_uniform_endpoint_remainder_221.json").write_text(
    json.dumps(out,indent=2),encoding="utf-8")
