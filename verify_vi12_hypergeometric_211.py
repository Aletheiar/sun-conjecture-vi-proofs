"""Exact rational pullback to a classical 3F2 equation, with target weights."""
import json
from pathlib import Path
import sympy as s
from sympy.functions.combinatorial.numbers import stirling
x,z,t=s.symbols("x z t")
root=Path(__file__).parent
data=json.loads((root/"vi12_rank19_operator_206.json").read_text())
L=s.sympify(data["euler"],locals={"x":x,"th":t})
zz=108*x*x*(1-4*x)
# W satisfies theta_z^3 W-z(theta_z+1/2)(theta_z+1/3)(theta_z+2/3)W=0.
Lw=s.expand(t**3-z*(t+s.Rational(1,2))*(t+s.Rational(1,3))*(t+s.Rational(2,3)))
awc=[s.factor(z**j*sum(Lw.coeff(t,k)*stirling(k,j,kind=2)
                       for k in range(j,4))) for j in range(4)]
C=s.zeros(3); C[0,1]=1; C[1,2]=1
for j in range(3): C[2,j]=s.factor(-awc[j]/awc[3])
# y=P(x) W(z(x)); row derivatives divided by P have rational coefficients.
plog=-2/(1-4*x)+3/(1-3*x)
N=s.diff(zz,x)*C.subs(z,zz)
rows=[s.Matrix([[1,0,0]])]
for j in range(3):
    rows.append((rows[-1].diff(x)+plog*rows[-1]+rows[-1]*N).applyfunc(s.cancel))
res=sum((x**j*sum(L.coeff(t,k)*stirling(k,j,kind=2)
                   for k in range(j,4))*rows[j] for j in range(4)),s.zeros(1,3))
assert all(s.factor(v)==0 for v in res)
assert s.factor(L.coeff(x,0)-t**3)==0
print("E1: exact pullback L(P W(z(x)))=0 and indicial theta^3.",flush=True)
# Unique analytic local solution, leading coefficient at index n is n^3.
# P(0)=W(0)=1. Therefore G=P W(z(x)), by E1 and the separately certified ODE.
normal=json.loads((root/"vi12_weight_translation_210.json").read_text())
out={}
for name,x0,ca,cb in [("VI1",s.Rational(1,36),17,66),("VI2",-s.Rational(1,60),31,126)]:
    A=9*x*x*(4*x-1)/((3*x-1)**2*(12*x-1))
    B=(4*x-1)*(6*x-1)/((3*x-1)*(12*x-1))
    aa=s.factor(ca+cb*A); bb=s.factor(cb*B)
    Pc=s.sqrt(1-4*x)/(1-3*x)
    cW=s.simplify((Pc*(aa+bb*x*plog)).subs(x,x0))
    cDW=s.simplify((Pc*bb*x*s.diff(zz,x)/zz).subs(x,x0))
    value=s.factor(zz.subs(x,x0))
    if name=="VI1":
        factor=80*s.sqrt(2)/11
        assert s.simplify(cW-2*factor)==0 and s.simplify(cDW-15*factor)==0
        ram=27/(4*s.pi); expected=540*s.sqrt(2)/(11*s.pi)
    else:
        factor=352*s.sqrt(15)/189
        assert s.simplify(cW-4*factor)==0 and s.simplify(cDW-33*factor)==0
        ram=15*s.sqrt(3)/(2*s.pi); expected=880*s.sqrt(5)/(21*s.pi)
    assert s.simplify(factor*ram-expected)==0
    out[name]={"z":str(value),"factor":str(factor),
               "classical_series_value":str(ram),"original_value":str(expected)}
    print(name,": z=",value,"factor=",factor,"original value=",expected,flush=True)
(root/"vi12_hypergeometric_translation_211.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
print("Uses the two classical Ramanujan 1/pi evaluations as external theorems.")
