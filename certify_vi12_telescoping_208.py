"""Try an exact telescoping certificate; failed ansatz is not a disproof."""
import json
from pathlib import Path
import sympy as s
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--target",choices=["normal","ode"],default="normal")
args=parser.parse_args()

x,n=s.symbols("x n")
b=1-6*x; d=-3*(1-4*x); w=x/(1-4*x); a=-3/(1-3*x)
M=s.zeros(4)
for j in range(4):
    M[j,j]=s.factor((3-j)*n*a+j*a*(n-1+(2*n-1)*b/d)+n/(x*(1-4*x)))
    if j<3: M[j,j+1]=(3-j)*n*a
    if j: M[j,j-1]=-j*n*a/d
N=s.zeros(4)
al=(2*n+1)*b/(n+1); be=-n*d/(n+1)
for j in range(4):
    for k in range(4-j):
        N[j,k]=s.binomial(3-j,k)*al**(3-j-k)*be**k*w

def simp(m):
    return m.applyfunc(s.cancel)

v0=s.Matrix([[1,0,0,0]])
v1=simp(v0*M)
pol=[-9*x*x*(4*x-1),-(3*x-1)*(4*x-1)*(6*x-1),0,(3*x-1)**2*(12*x-1)]
E=simp(pol[0]*v0+pol[1]*x*v1+pol[3]*n*v0)
if args.target=="ode":
    from sympy.functions.combinatorial.numbers import stirling
    data=json.loads(Path(__file__).with_name("vi12_rank19_operator_206.json").read_text())
    th=s.symbols("th")
    L=s.sympify(data["euler"],locals={"x":x,"th":th})
    rows=[v0,v1]
    rows.append(simp(rows[-1].diff(x)+rows[-1]*M))
    rows.append(simp(rows[-1].diff(x)+rows[-1]*M))
    E=simp(sum(((x**j*sum(L.coeff(th,k)*stirling(k,j,kind=2)
                         for k in range(j,4)))*rows[j]
               for j in range(4)),s.zeros(1,4)))
print(args.target,"source ready.",flush=True)
# Solve over Q(x) by polynomial identities in n. First eliminate rational n poles.
found_cert=False
for degree in range(1,5):
    cs=s.symbols("c:"+str(4*(degree+1)))
    R=s.Matrix([[sum(cs[j*(degree+1)+i]*n**i for i in range(degree+1))
                 for j in range(4)]])
    residual=simp(R.subs(n,n+1)*N-R-E)
    eq=[]
    for entry in residual:
        eq.extend(s.Poly(s.together(entry).as_numer_denom()[0],n).all_coeffs())
    print("degree",degree,"unknowns",len(cs),"equations",len(eq),flush=True)
    mat,rhs=s.linear_eq_to_matrix(eq,cs)
    # Solve in the rational function field with fraction-free DomainMatrix.
    from sympy.polys.matrices import DomainMatrix
    from sympy.polys.domains import QQ
    field=QQ.frac_field(x)
    dm=DomainMatrix.from_Matrix(mat.row_join(rhs)).convert_to(field)
    rr,piv=dm.rref()
    rm=rr.to_Matrix()
    if len(cs) in piv:
        print("No polynomial certificate of this degree.",flush=True)
        continue
    sol={c:s.Integer(0) for c in cs}
    for i,j in enumerate(piv):
        sol[cs[j]]=rm[i,-1]
    cert=simp(R.subs(sol))
    assert all(v==0 for v in simp(cert.subs(n,n+1)*N-cert-E))
    print("Exact certificate:",list(cert),flush=True)
    print("n=0 boundary:",list(cert.subs(n,0)),flush=True)
    Path(__file__).with_name("vi12_telescoping_208.json").write_text(
        json.dumps({"relation":args.target,"certificate":[str(c) for c in cert]},
                   indent=2),encoding="utf-8")
    found_cert=True
    break
if not found_cert:
    raise SystemExit("No certificate found in the searched polynomial ansatz.")
