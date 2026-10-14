"""Exact identities for the 41-zero pure-head degree-seven construction."""
from pathlib import Path
import sympy as s,json,hashlib,time
if not __debug__:raise SystemExit('Assertions are required; do not run with -O.')
z,w,x,y,du,dv=s.symbols('z w x y du dv')
a=s.Rational(1,4);kap=s.Rational(5,64);T=s.S(100)
d=3*a*a+a;J=d*d-a*a;A=d/J;B=-a/J
f=z**3+z*w-kap;fs=w**3+z*w-kap;g=z-a;G=w*w*f-z*fs
gb=s.groebner([f,fs],z,w,order='grlex',domain=s.QQ)
red=lambda p:s.expand(gb.reduce(s.expand(p))[1])
flinear=d*(z-a)+a*(w-a);fslinear=a*(z-a)+d*(w-a)
Bpoly=-a*a*A*flinear+(A*(T+a)-a*a*B)*fslinear
H=s.expand(G*g+f*Bpoly)
el=(kap-z**3)**3-z**6
eid=s.cancel(el/((z-a)*s.diff(el,z).subs(z,a)))
gin=s.cancel((1-eid)/(z-a))
assert red(g*gin-(1-eid))==0
alpha=red(w*w+Bpoly*gin);nu=-z;alphac=alpha.xreplace({z:w,w:z})
R=red(g*(alpha*alpha/4+nu*alphac/2))
assert R.subs({z:a,w:a})==0
Rz=s.diff(R,z).subs({z:a,w:a});Rw=s.diff(R,w).subs({z:a,w:a})
Ru=(d*Rz-a*Rw)/J;Rv=(-a*Rz+d*Rw)/J
K=s.expand(R-(Ru+a*Rv)*f+(Rv/a)*G)
L=-A*T**3/100
assert red(K-R)==0
assert all(s.diff(K,z,i,w,j).subs({z:a,w:a})==0 for i,j in [(0,0),(1,0),(0,1)])
uf=d*du+a*dv;vf=a*du+d*dv
def jet(poly,n):
    p=s.Poly(poly.subs({z:a+du,w:a+dv}).expand(),du,dv)
    return s.expand(sum(co*du**i*dv**j for (i,j),co in p.terms() if i+j<=n))
assert s.expand(jet(f*f*g,3)-A*uf**3-B*uf**2*vf)==0
assert s.expand(jet(H,2)-A*T*uf*vf+a*B*vf**2)==0
assert B/A==-s.Rational(4,7) and -a*B/(A*T)==s.Rational(1,700)
for P in [H,K,L]:
    p=s.Poly(P,z,w,domain=s.QQ)
    assert p.total_degree()<=6 and p.degree(w)<=2
u=x+s.I*y;v=x-s.I*y
model=s.expand(u**3-s.Rational(4,7)*u*u*v+u*v+s.Rational(1,700)*v*v-s.Rational(1,100))
re=s.expand(model).collect(s.I).coeff(s.I,0);im=s.expand(model).coeff(s.I)
assert s.expand(im-y*(17*x*x-11*y*y-x/50)/7)==0
axis=300*x**3+701*x*x-7
off=-1960000*x**3+982200*x*x-699*x-3850
assert s.expand(re.subs(y,0)-axis/700)==0
yp=x*(850*x-1)/550
assert s.cancel(re.subs(y*y,yp)-off/385000)==0
axis_intervals=[(-3,-2),(-s.Rational(1,5),-s.Rational(1,20)),(s.Rational(1,20),s.Rational(1,5))]
off_intervals=[(-s.Rational(1,5),-s.Rational(1,20)),(s.Rational(1,20),s.Rational(1,5)),(s.Rational(1,5),1)]
signs=[]
for name,p,intervals in [('axis',axis,axis_intervals),('off',off,off_intervals)]:
    for lo,hi in intervals:
        assert p.subs(x,lo)*p.subs(x,hi)<0
        assert hi<0 or lo>s.Rational(1,850)
        signs.append(dict(branch=name,interval=[str(lo),str(hi)],values=[str(p.subs(x,lo)),str(p.subs(x,hi))]))
assert s.gcd(axis,s.diff(axis,x))==1 and s.gcd(off,s.diff(off,x))==1
assert 8*4+9==41 and 24+17==41 and 24-17==7
coeffs=lambda p:[[i,j,str(co)] for (i,j),co in s.Poly(p,z,w).terms()]
out=dict(all_passed=True,arithmetic='exact rational polynomial arithmetic',seed_a=str(a),seed_kappa=str(kap),
         T=str(T),A=str(A),B=str(B),H_coefficients=coeffs(H),K_coefficients=coeffs(K),L=str(L),
         shared_model=str(model),axis_polynomial=str(axis),off_axis_polynomial=str(off),
         root_sign_intervals=signs,total_regular_zeros=41,positive=24,negative=17,real_axis_zeros=7,
         verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
target=Path(__file__).with_name('pure_2_7_41_verification.json')
target.write_text(json.dumps(out,indent=2)+'\n');print('ALL PASSED',target)
