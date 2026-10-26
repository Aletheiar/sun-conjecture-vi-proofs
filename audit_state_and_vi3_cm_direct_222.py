"""Independent formal-state and direct K,E elimination audit.

No hypergeometric Ramanujan R3 formula or final-translation JSON is read.
Classical CM(7) modulus, period ratio and alpha(7) are explicit inputs.
"""
from pathlib import Path
import json
import sympy as s
n,b,d,w,sigma,tau,omega=s.symbols("n b d w sigma tau omega")
U,V=s.symbols("U V")
basis=[U**(3-j)*V**j for j in range(4)]
Up=n*(sigma*U+tau*V)
Vp=-n*tau/d*U+((n-1)*sigma+(2*n-1)*b*tau/d)*V
M=s.zeros(4); N=s.zeros(4)
for j,Bj in enumerate(basis):
    chain=s.Poly(s.expand(s.diff(Bj,U)*Up+s.diff(Bj,V)*Vp+n*omega*Bj),U,V)
    for k,Bk in enumerate(basis):
        M[j,k]=s.factor(chain.coeff_monomial(Bk))
        if abs(j-k)>1: assert M[j,k]==0
    assert s.factor(M[j,j]-((3-j)*n*sigma+
          j*((n-1)*sigma+(2*n-1)*b*tau/d)+n*omega))==0
    if j<3: assert M[j,j+1]==(3-j)*n*tau
    if j: assert s.factor(M[j,j-1]+j*n*tau/d)==0
    shifted=s.Poly(s.expand(w*((2*n+1)*b*U/(n+1)-n*d*V/(n+1))**(3-j)*U**j),U,V)
    for k,Bk in enumerate(basis):
        N[j,k]=s.factor(shifted.coeff_monomial(Bk))
        predicted=(w*s.binomial(3-j,k)*((2*n+1)*b/(n+1))**(3-j-k)
                   *(-n*d/(n+1))**k) if k<=3-j else 0
        assert s.factor(N[j,k]-predicted)==0
print("E1: all-index state matrices derived from formal chain rule and binomial expansion.",flush=True)
# No R3 substitution: directly eliminate the elliptic values K,E,Kc,Ec.
t=s.symbols("t",positive=True)
K,E,Kc,Ec=s.symbols("K E Kc Ec",nonzero=True)
mm=(1-s.sqrt(t))**2/(2*(1+t))
md=s.diff(mm,t)
P=s.sqrt(1-t)/(1+t-t*t)
AA=3*(t-1)*(t*t+1)*(t*t+t+1)/(8*(t*t-t-1)**2)
BB=(t-1)*(t+1)*(t*t+2*t-1)/(4*(t*t-t-1))
RR=(t*t+1)*(t+1)*s.sqrt(1-t)/(s.sqrt(t)*(t*t-t-1)**2)
G=P*(3*K*K+Kc*Kc)/s.pi**2
Kdot=md*(E/(2*mm*(1-mm))-K/(2*mm))
Kcdot=-md*(Ec/(2*mm*(1-mm))-Kc/(2*(1-mm)))
Gdot=s.diff(P,t)*(3*K*K+Kc*Kc)/s.pi**2+2*P*(3*K*Kdot+Kc*Kcdot)/s.pi**2
t0=s.Rational(7,9)
assert s.simplify(mm.subs(t,t0)-(8-3*s.sqrt(7))/16)==0
alpha7=(s.sqrt(7)-2)/2
# alpha(7)=pi/(4K^2)+sqrt(7)*(1-E/K).
E7=(1-alpha7/s.sqrt(7))*K+s.pi/(4*s.sqrt(7)*K)
Kc7=s.sqrt(7)*K
Ec7=(s.pi/2+K*Kc7-E7*Kc7)/K
assert s.simplify(K*Ec7+E7*Kc7-K*Kc7-s.pi/2)==0
assert s.simplify(s.pi/(4*K*K)+s.sqrt(7)*(1-E7/K)-alpha7)==0
elliptic_sub={E:E7,Kc:Kc7,Ec:Ec7}
lhs=((1147+3990*AA)*G+3990*BB*Gdot+3990*3*RR/(4*s.pi)).subs(t,t0)
pre_alpha=s.expand(lhs.subs({Kc:Kc7,Ec:(s.pi/2+K*Kc7-E*Kc7)/K},simultaneous=True))
coef_K2=s.simplify(pre_alpha.coeff(K,2))
coef_KE=s.simplify(pre_alpha.coeff(K,1).coeff(E,1))
constant=s.simplify(pre_alpha-coef_K2*K*K-coef_KE*K*E)
assert not ({K,E}&constant.free_symbols)
print("Direct coefficients [K^2,K*E,constant]:",[coef_K2,coef_KE,constant],flush=True)
value=s.simplify(s.expand(lhs.subs(elliptic_sub,simultaneous=True)))
target=432*(195*s.sqrt(14)+94*s.sqrt(2))/(95*s.pi)
assert s.simplify(value-target)==0
assert K not in value.free_symbols
print("E2: direct elliptic CM(7) elimination leaves no K-dependent term.",flush=True)
print("E3: final value equals",s.simplify(value),flush=True)
missing=s.simplify((lhs-3990*3*RR.subs(t,t0)/(4*s.pi)).subs(elliptic_sub,simultaneous=True)-target)
assert missing!=0
assert s.simplify(missing+3990*3*RR.subs(t,t0)/(4*s.pi))==0
print("E4: omitting the nonzero boundary term fails the exact target.",flush=True)
Path(__file__).with_name("vi3_direct_cm_audit_222.json").write_text(
    json.dumps({"method":"formal generic state derivation; direct K E Kc Ec elimination",
                "classical_inputs":["m7=(8-3sqrt7)/16","Kc/K=sqrt7","alpha7=(sqrt7-2)/2"],
                "R3_used":False,"target":str(target),
                "K2_coefficient":str(coef_K2),"KE_coefficient":str(coef_KE),
                "constant_after_legendre":str(constant),
                "omitted_boundary_error":str(missing)},indent=2),encoding="utf-8")
