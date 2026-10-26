"""Independent exact checks for the VI3 proof; numerical checks are optional.

Analytic domination and branch matching are proved in the accompanying Markdown.
Classical CM(7) data are explicit external inputs, not numerical discoveries.
"""
import json
from pathlib import Path
import sympy as s
t=s.symbols("t",positive=True)
n,v=s.symbols("n v")
root=Path(__file__).parent
def parse(st): return s.sympify(st,locals={"t":t,"n":n})
def simp(m): return m.applyfunc(s.cancel)
def zero(m): assert all(s.factor(e)==0 for e in m)
b=2*(t*t+t-1); c=(1+t-t*t)**2; d=s.factor(b*b-4*c)
w=1/(64*(t-1)*(t+1)**2)
sigma=s.factor((b*s.diff(b,t)-s.diff(d,t)/2)/(b*b-d))
tau=s.factor(-(d*s.diff(b,t)-b*s.diff(d,t)/2)/(b*b-d))
den=1-2*b*v+d*v*v
assert s.factor((s.diff(b,t)*v-s.diff(d,t)*v*v/2)/den
     -sigma*v*(b-d*v)/den-tau*v*(v*(b-d*v)/den+1))==0
M=s.zeros(4); N=s.zeros(4)
for j in range(4):
    M[j,j]=(3-j)*n*sigma+j*((n-1)*sigma+(2*n-1)*b*tau/d)+n*s.diff(w,t)/w
    if j<3: M[j,j+1]=(3-j)*n*tau
    if j: M[j,j-1]=-j*n*tau/d
    for k in range(4-j):
        N[j,k]=w*s.binomial(3-j,k)*((2*n+1)*b/(n+1))**(3-j-k)*(-n*d/(n+1))**k
M=simp(M); N=simp(N)
rows=[s.Matrix([[1,0,0,0]])]
for j in range(3): rows.append(simp(rows[-1].diff(t)+rows[-1]*M))
data=json.loads((root/"vi3_locus_operator_217.json").read_text())
op=[parse(a) for a in data["operator"]]
R=s.Matrix([[parse(a) for a in data["certificate"]]])
E=sum((op[j]*rows[j] for j in range(4)),s.zeros(1,4))
zero(R.subs(n,n+1)*N-R-E); zero(R.subs(n,0))
assert op[3]==1
print("E1: all-n generating-function identity, third-order certificate, zero boundary.",flush=True)
normal=json.loads((root/"vi3_normal_system_219.json").read_text())[0]
co=[parse(a) for a in normal["relation"]]
R2=s.Matrix([[parse(a) for a in normal["certificate"]]])
terms=rows[:3]+[n*rows[0],n*rows[1]]
E2=sum((co[j]*terms[j] for j in range(5)),s.zeros(1,4))
zero(R2.subs(n,n+1)*N-R2-E2); zero(R2.subs(n,0))
assert co[4]==1
AA=3*(t-1)*(t*t+1)*(t*t+t+1)/(8*(t*t-t-1)**2)
BB=(t-1)*(t+1)*(t*t+2*t-1)/(4*(t*t-t-1))
assert s.factor(BB+co[2])==0
assert s.factor(AA+s.diff(BB,t)+co[3]*BB+co[1])==0
assert s.factor(s.diff(AA,t)+co[3]*AA+co[0])==0
logRR=2*t/(t*t+1)+1/(t+1)+1/(2*(t-1))-1/(2*t)-2*(2*t-1)/(t*t-t-1)
assert s.factor(logRR+co[3])==0
assert AA.subs(t,0)==-s.Rational(3,8)
assert BB.subs(t,0)==-s.Rational(1,4)
print("E2: normal derivative system, rational particular solution and homogeneous factor.",flush=True)
p=s.factor(op[2]/3)
q=s.factor((op[1]-s.diff(p,t)-2*p*p)/4)
assert s.factor(op[0]-2*s.diff(q,t)-4*p*q)==0
mm=(1-s.sqrt(t))**2/(2*(1+t))
zz=(1-t)**2/(1+t)**2
assert s.simplify(zz-4*mm*(1-mm))==0
lp=1/(2*(t-1))-(2*t-1)/(t*t-t-1)
md=s.diff(mm,t)
pk=s.simplify((1-2*mm)/(mm*(1-mm))*md-s.diff(md,t)/md)
qk=s.simplify(-md*md/(4*mm*(1-mm)))
assert s.factor(p+lp-pk)==0
assert s.factor(q+p*lp/2+s.diff(lp,t)/2+lp*lp/4-qk)==0
print("E3: exact Legendre pullback; P*K(m)^2, P*K(m)K(1-m), P*K(1-m)^2 solve L.",flush=True)
# Asymptotic normalization constants, using Gauss summation and Legendre at m=1/2.
gamma=s.symbols("gamma",positive=True)
K0=gamma**2/(4*s.sqrt(s.pi))
K1=s.pi/(2*K0)
h0=4*K0*K0/s.pi**2
c0=16*s.pi/gamma**4
assert s.simplify(c0-4*K1*K1/s.pi**2)==0
assert s.simplify(h0*c0-4/s.pi**2)==0
# C_n^minus=-C_n^plus/(2n-1), C_n^contiguous=(2n+1)*C_n^minus.
an=s.symbols("an")
bn=-an/(2*n-1)
assert s.factor((2*n+1)*bn-2*bn+an)==0
CH=(K1*K1+K0*K0)/s.pi**2
source_linear=CH-3*(h0-c0)/4
candidate_linear=4*CH-3*h0/2
assert s.simplify(source_linear-candidate_linear)==0
assert s.simplify(source_linear-(-h0/2+c0))==0
assert s.simplify(-4*K0*K1/s.pi**2+2/s.pi)==0
# H ~ 1/(pi sqrt(t)), G' ~ -1/(pi sqrt(t)).
assert 1-BB.subs(t,0)*(-1)==s.Rational(3,4)
print("E4: exact coefficients of 1,sqrt(t),t; normal constant 3/(4*pi).",flush=True)
t0=s.Rational(7,9)
P=s.sqrt(1-t)/(1+t-t*t)
RR=(t*t+1)*(t+1)*s.sqrt(1-t)/(s.sqrt(t)*(t*t-t-1)**2)
assert s.simplify(mm.subs(t,t0)-(8-3*s.sqrt(7))/16)==0
assert zz.subs(t,t0)==s.Rational(1,64)
assert s.simplify(P.subs(t,t0)-27*s.sqrt(2)/95)==0
a=1147+3990*AA; bb=3990*BB
# CM ratio K(1-m0)/K(m0)=sqrt(7), plus Legendre Wronskian.
cH=s.simplify((s.Rational(5,2)*P*(a+bb*lp)).subs(t,t0))
cD=s.simplify((s.Rational(5,2)*P*bb*s.diff(zz,t)/zz).subs(t,t0))
extra=s.simplify((-bb*P*s.sqrt(7)*md/(2*mm*(1-mm))
                  +3990*s.Rational(3,4)*RR).subs(t,t0))
assert s.simplify(cH/5-cD/42)==0
target=s.Rational(432,95)*(195*s.sqrt(14)+94*s.sqrt(2))
assert s.simplify(16*cH/5+extra-target)==0
print("E5: CM(7) reduction and final target algebra: exact.",flush=True)
print("Classical input: CM(7) ratio and (5+42 theta)3F2(1/64)=16/pi.",flush=True)
out={"G_formula":"P*(3*K(m)^2+K(1-m)^2)/pi^2",
     "normal_A":str(s.factor(AA)),"normal_B":str(s.factor(BB)),
     "normal_R":str(RR),"normal_constant":"3/(4*pi)",
     "classical_H_coefficient":str(cH),"classical_thetaH_coefficient":str(cD),
     "extra_pi_coefficient":str(extra),"target_pi_coefficient":str(target)}
(root/"vi3_verified_formula_220.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
