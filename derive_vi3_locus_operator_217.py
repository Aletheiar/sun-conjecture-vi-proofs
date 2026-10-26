"""Seek an exact order-three telescoper on the new VI3 family."""
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ
from pathlib import Path
import json
t,n=s.symbols("t n")
b=2*(t*t+t-1); c=(1+t-t*t)**2
d=s.factor(b*b-4*c)
w=1/(64*(t-1)*(t+1)**2)
sigma=s.factor((b*s.diff(b,t)-s.diff(d,t)/2)/(b*b-d))
tau=s.factor(-(d*s.diff(b,t)-b*s.diff(d,t)/2)/(b*b-d))
M=s.zeros(4); N=s.zeros(4)
for j in range(4):
    M[j,j]=s.factor((3-j)*n*sigma+j*((n-1)*sigma+(2*n-1)*b*tau/d)+n*s.diff(w,t)/w)
    if j<3: M[j,j+1]=(3-j)*n*tau
    if j: M[j,j-1]=-j*n*tau/d
    for k in range(4-j):
        N[j,k]=w*s.binomial(3-j,k)*((2*n+1)*b/(n+1))**(3-j-k)*(-n*d/(n+1))**k
def simp(m): return m.applyfunc(s.cancel)
rows=[s.Matrix([[1,0,0,0]])]
for j in range(3):
    rows.append(simp(rows[-1].diff(t)+rows[-1]*M))
print("Derivative sources built.",flush=True)
degree=3; order=3
cs=s.symbols("c:"+str(4*(degree+1)))
ls=s.symbols("l:"+str(order+1))
R=s.Matrix([[sum(cs[j*(degree+1)+i]*n**i for i in range(degree+1)) for j in range(4)]])
target=sum((ls[j]*rows[j] for j in range(order+1)),s.zeros(1,4))
eqmat=simp(R.subs(n,n+1)*N-R-target)
eq=[]
for ent in eqmat:
    eq.extend(s.Poly(s.together(ent).as_numer_denom()[0],n).all_coeffs())
mat,_=s.linear_eq_to_matrix(eq,cs+ls)
print("Matrix",mat.shape,flush=True)
dm=DomainMatrix.from_Matrix(mat).convert_to(QQ.frac_field(t))
ker=dm.nullspace().to_Matrix()
print("Nullity",ker.rows,flush=True)
for i in range(ker.rows):
    vec=[s.cancel(v) for v in ker.row(i)]
    if vec[-1]==0: continue
    norm=vec[-1]; vec=[s.cancel(v/norm) for v in vec]
    sol=dict(zip(cs+ls,vec))
    cert=simp(R.subs(sol))
    op=[s.factor(v) for v in vec[len(cs):]]
    assert all(v==0 for v in simp(eqmat.subs(sol)))
    assert all(v==0 for v in cert.subs(n,0))
    print("Operator:",op,flush=True)
    Path(__file__).with_name("vi3_locus_operator_217.json").write_text(
        json.dumps({"operator":[str(v) for v in op],
                    "certificate":[str(v) for v in cert],
                    "sigma":str(sigma),"tau":str(tau)},indent=2),encoding="utf-8")
    break
else:
    raise SystemExit("No order-three polynomial certificate in this ansatz.")
