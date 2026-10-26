"""Exact joint period relation for G,G',G'',H,H' on the new VI3 family."""
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ
from pathlib import Path
import json
t,n=s.symbols("t n")
b=2*(t*t+t-1); c=(1+t-t*t)**2; d=s.factor(b*b-4*c)
w=1/(64*(t-1)*(t+1)**2)
sigma=s.factor((b*s.diff(b,t)-s.diff(d,t)/2)/(b*b-d))
tau=s.factor(-(d*s.diff(b,t)-b*s.diff(d,t)/2)/(b*b-d))
M=s.zeros(4); N=s.zeros(4)
for j in range(4):
    M[j,j]=(3-j)*n*sigma+j*((n-1)*sigma+(2*n-1)*b*tau/d)+n*s.diff(w,t)/w
    if j<3: M[j,j+1]=(3-j)*n*tau
    if j: M[j,j-1]=-j*n*tau/d
    for k in range(4-j):
        N[j,k]=w*s.binomial(3-j,k)*((2*n+1)*b/(n+1))**(3-j-k)*(-n*d/(n+1))**k
def simp(m): return m.applyfunc(s.cancel)
e0=s.Matrix([[1,0,0,0]])
e1=simp(e0*M); e2=simp(e1.diff(t)+e1*M)
terms=[e0,e1,e2,n*e0,n*e1]
degree=2
cs=s.symbols("c:"+str(4*(degree+1))); ls=s.symbols("l:5")
R=s.Matrix([[sum(cs[j*(degree+1)+i]*n**i for i in range(degree+1)) for j in range(4)]])
target=sum((ls[j]*terms[j] for j in range(5)),s.zeros(1,4))
eqmat=simp(R.subs(n,n+1)*N-R-target)
eq=[]
for ent in eqmat: eq.extend(s.Poly(s.together(ent).as_numer_denom()[0],n).all_coeffs())
mat,_=s.linear_eq_to_matrix(eq,cs+ls)
print("Joint matrix",mat.shape,flush=True)
ker=DomainMatrix.from_Matrix(mat).convert_to(QQ.frac_field(t)).nullspace().to_Matrix()
print("Joint nullity",ker.rows,flush=True)
out=[]
for i in range(ker.rows):
    vec=[s.cancel(v) for v in ker.row(i)]
    if not any(vec[len(cs):]): continue
    norm=next(v for v in reversed(vec[len(cs):]) if v!=0)
    vec=[s.cancel(v/norm) for v in vec]
    sol=dict(zip(cs+ls,vec))
    cert=simp(R.subs(sol))
    op=[s.factor(v) for v in vec[len(cs):]]
    assert all(v==0 for v in simp(eqmat.subs(sol)))
    assert all(v==0 for v in cert.subs(n,0))
    print("Relation G,G',G'',H,H':",op,flush=True)
    out.append({"relation":[str(v) for v in op],"certificate":[str(v) for v in cert]})
Path(__file__).with_name("vi3_normal_system_219.json").write_text(
    json.dumps(out,indent=2),encoding="utf-8")
if not out: raise SystemExit("No joint relation in this ansatz.")
