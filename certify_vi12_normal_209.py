"""Find a telescoper for the residual of the normal-derivative identity."""
import sys
sys.argv=[sys.argv[0]]
# Share exactly checked source matrices without re-solving polynomial certificates.
from pathlib import Path
ns={}
source=Path(__file__).with_name("certify_vi12_telescoping_208.py").read_text()
exec(source[:source.index("# Solve over Q(x)")],ns)
s,x,n,M,N,E,simp=[ns[k] for k in ["s","x","n","M","N","E","simp"]]
from sympy.polys.matrices import DomainMatrix
from sympy.polys.domains import QQ
import json
terms=[E]
for order in range(1,5):
    terms.append(simp(terms[-1].diff(x)+terms[-1]*M))
    degree=order+1
    cs=s.symbols("c:"+str(4*(degree+1)))
    ls=s.symbols("l:"+str(order+1))
    R=s.Matrix([[sum(cs[j*(degree+1)+i]*n**i for i in range(degree+1))
                 for j in range(4)]])
    target=sum((ls[j]*terms[j] for j in range(order+1)),s.zeros(1,4))
    eqmat=simp(R.subs(n,n+1)*N-R-target)
    eq=[]
    for ent in eqmat:
        eq.extend(s.Poly(s.together(ent).as_numer_denom()[0],n).all_coeffs())
    mat,_=s.linear_eq_to_matrix(eq,cs+ls)
    print("order",order,"matrix",mat.shape,flush=True)
    dm=DomainMatrix.from_Matrix(mat).convert_to(QQ.frac_field(x))
    ker=dm.nullspace().to_Matrix()
    found=False
    for i in range(ker.rows):
        vec=[s.cancel(v) for v in ker.row(i)]
        if not any(vec[len(cs):]): continue
        norm=vec[-1]
        if norm==0: norm=next(v for v in vec[len(cs):] if v!=0)
        vec=[s.cancel(v/norm) for v in vec]
        sol=dict(zip(cs+ls,vec))
        cert=simp(R.subs(sol))
        op=[s.factor(v) for v in vec[len(cs):]]
        assert all(v==0 for v in simp(eqmat.subs(sol)))
        print("Residual ODE coefficients:",op,flush=True)
        print("Boundary:",list(cert.subs(n,0)),flush=True)
        Path(__file__).with_name("vi12_normal_residual_209.json").write_text(
            json.dumps({"operator":[str(v) for v in op],
                        "certificate":[str(v) for v in cert]},indent=2),encoding="utf-8")
        found=True
        break
    if found: break
else:
    raise SystemExit("No residual differential certificate found.")
