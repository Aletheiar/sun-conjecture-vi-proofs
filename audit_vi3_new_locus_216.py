"""A new rational parameter locus passing through VI3, with actual section."""
import json
from pathlib import Path
import sympy as s
t,u,l,k=s.symbols("t u l k")
lam=t-1/t
kk=t**3/((t-1)*(t+1)**2)
A=s.factor((1-kk)/(1-lam**2*kk))
B=s.factor((1-lam*kk)/(1-lam**2*kk))
D=s.factor((1-lam**3*kk)/(1-lam**2*kk))
assert s.factor(A*(2*B*D-1)+B**2-B**3*D)==0
assert s.factor(kk**2*lam**3-kk*lam**3+kk*lam**2-2*kk*lam+1)==0
T=-(u+D)/D
alpha=(B*u+1)**2
beta=(A*u+B)*(u+D)
V2=-A/(D*B**2)*u*(u+D)**2*(B**2*u+2*B-1/D)**2
assert s.factor(T*(T+alpha)*(T+beta)-V2)==0
b=2*(t*t+t-1); c=(1+t-t*t)**2
w=1/(64*(t-1)*(t+1)**2)
assert s.factor((b+2*(1+t-t*t))**3*w-kk)==0
assert s.factor((b-2*(1+t-t*t))/(b+2*(1+t-t*t))-lam)==0
t0=s.Rational(7,9); scale=s.Integer(81)
assert s.factor((scale*b).subs(t,t0)-62)==0
assert s.factor((scale**2*c).subs(t,t0)-9025)==0
assert s.factor((w/scale**3).subs(t,t0)+s.Rational(1,23887872))==0
assert A.subs(t,t0)==s.Rational(729,512)
assert B.subs(t,t0)==s.Rational(9,16)
assert D.subs(t,t0)==s.Rational(7,9)
print("New family:",b,c,w,flush=True)
print("A B D:",A,B,D,flush=True)
print("Exact section identity and VI3 scaling verified.",flush=True)
print("No independence, rank or period special value claimed.",flush=True)
Path(__file__).with_name("vi3_new_locus_216.json").write_text(
    json.dumps({"b":str(b),"c":str(c),"w":str(w),"A":str(A),"B":str(B),"D":str(D),
                "parameter":"t","VI3_t":"7/9","VI3_scale":"81"},indent=2),
    encoding="utf-8")
