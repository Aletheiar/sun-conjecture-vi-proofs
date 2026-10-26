"""Exact symmetric-square and bounded Gauss pullback recognition."""
import json
from pathlib import Path
import sympy as s
t,z=s.symbols("t z")
root=Path(__file__).parent
d=json.loads((root/"vi3_locus_operator_217.json").read_text())
c=[s.sympify(v) for v in d["operator"]]
p=s.factor(c[2]/3)
q=s.factor((c[1]-s.diff(p,t)-2*p*p)/4)
assert s.factor(c[0]-2*s.diff(q,t)-4*p*q)==0
print("Exact symmetric square verified.",flush=True)
candidates=[((1-t)/(1+t))**2,4*t/(1+t)**2,1-t*t,
            4*t*(1-t),(1-t)/(1+t),t*t,(1-t)**2,
            -4*t/(1-t)**2,4*t*t/(1+t*t)**2]
found=[]
for zz in candidates:
    dz=s.diff(zz,t)
    pp=s.factor(((1-s.Rational(3,2)*z)/(z*(1-z))).subs(z,zz)*dz-s.diff(dz,t)/dz)
    ll=s.factor((pp-p)/2)
    for aa in [s.Rational(1,2),s.Rational(1,3),s.Rational(1,4),s.Rational(1,6)]:
        qq=s.factor((-aa*(1-aa)/(4*z*(1-z))).subs(z,zz)*dz*dz)
        obs=s.factor(q+p*ll+s.diff(ll,t)+ll*ll-qq)
        if obs==0:
            entry={"argument":str(s.factor(zz)),"a":str(aa),
                   "log_h_gauge":str(ll),"log_G_gauge":str(s.factor(2*ll)),
                   "G_gauge_integral":str(s.integrate(2*ll,t))}
            print("MATCH:",entry,flush=True)
            found.append(entry)
        else:
            print("No match:",str(s.factor(zz)),str(aa),flush=True)
(root/"vi3_hypergeometric_218.json").write_text(
    json.dumps({"p":str(p),"q":str(q),"matches":found},indent=2),encoding="utf-8")
