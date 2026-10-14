"""Exact algebraic checks for the unique deltoid incircle and quartic disks.
The global maximal-circle argument is analytic in the accompanying proof;
these exact identities check its algebra and the parameter transformations.
"""
from pathlib import Path
import sympy as s
import hashlib,json
if not __debug__:raise RuntimeError('Run without -O: proof assertions must remain enabled')
q,h,C,Cb,R,S,Sb,U,a,b,t,x=s.symbols('q h C Cb R S Sb U a b t x')
D=lambda z,zb:z*z*zb*zb+18*z*zb-4*(z**3+zb**3)-27
Delta=4*a**3+(36*a-4)*b*b-27*(a*a+b*b)**2
assert s.expand(D(3-27*(a+s.I*b),3-27*(a-s.I*b))+27**3*Delta)==0
assert s.expand(R**4+8*R**3+18*R**2-27-(R-1)*(R+3)**3)==0
assert s.expand(R**4-8*R**3+18*R**2-27-(R-3)**3*(R+1))==0
gamma=2*q+q**-2;gb=2/q+q*q
H=2*q**6-C*q**5-2*Cb*q**4+(5+C*Cb-R*R)*q**3-2*C*q**2-Cb*q+2
assert s.expand(q**3*((gamma-C)*(gb-Cb)-R*R)-H)==0
factor=2*(q**3-S*q**2+U*Sb*q-U)**2
poly=s.Poly(s.expand(factor-H.subs({C:4*S,Cb:4*Sb})),q)
assert poly.nth(5)==0
assert s.expand(poly.nth(4)-2*(S*S+(4+2*U)*Sb))==0
assert s.rem(poly.nth(1),U*U-1,U)==0
assert s.expand(factor.subs({S:0,Sb:0})-2*(q**3-U)**2)==0
for u in [s.S(-1),s.S(1)]:
    assert 4+2*u>0
    assert (4*(4+2*u)) in [8,24]
    assert (5+4*u) in [1,9]
derivative=2*s.I*(h*h-h**-4)
tangent=-4/h*(h**3-h**-3)/(2*s.I)
assert s.expand(derivative-tangent)==0
theta=s.symbols('theta',real=True)
series=s.series(2*s.exp(s.I*theta)+s.exp(-2*s.I*theta),theta,0,4).removeO()
assert s.expand(series-(3-3*theta**2+s.I*theta**3))==0
assert D(s.S(-2),s.S(-2))>0
psi=x**3-x*x+2*a*x-a*a-b*b
assert s.expand(s.discriminant(psi,x)-Delta)==0
at=t-s.Rational(3,2)*t*t;bt2=t**3*(1-s.Rational(9,4)*t)
radial=s.expand(psi.subs(a,at).subs(b*b,bt2))
assert s.expand(radial-(x-t)**2*(x-(1-2*t)))==0
assert s.expand(Delta.subs(a,at).subs(b*b,bt2))==0
distance=s.expand((at-a)**2+bt2)
expected=a*a-2*a*t+(1+3*a)*t*t-2*t**3
assert s.expand(distance-expected)==0
assert s.expand(s.diff(distance,t)+2*(3*t-1)*(t-a))==0
assert s.expand(distance.subs(t,a)-a**3)==0
assert s.expand(distance.subs(t,s.Rational(4,9))-(a-s.Rational(4,27))**2)==0
assert s.expand(h**3+h*h-s.Rational(4,27)-(h-s.Rational(1,3))*(h+s.Rational(2,3))**2)==0
contacts=[]
for c in [-s.S(1),(1+s.I*s.sqrt(3))/2,(1-s.I*s.sqrt(3))/2]:
    assert s.simplify(c**3+1)==0
    tau=s.simplify((3-c)/27)
    assert s.simplify((tau-s.Rational(1,9))*s.conjugate(tau-s.Rational(1,9))-s.Rational(1,729))==0
    assert s.simplify(Delta.subs({a:s.re(tau),b:s.im(tau)}))==0
    assert s.simplify(tau)!=0
    for cusp in [s.Rational(1,6)+s.I/(6*s.sqrt(3)),s.Rational(1,6)-s.I/(6*s.sqrt(3))]:assert s.simplify(tau-cusp)!=0
    contacts.append(str(tau))
out={'all_passed':True,'arithmetic':'exact symbolic identities',
     'deltoid_unique_maximal_disk':{'center':'0','radius':'1'},
     'folding_orders':'every integer k>=2; sharp all-center radius 1',
     'quartic_unique_maximal_disk':{'center':'1/9','radius':'1/27','exceptional_boundary_parameters':contacts},
     'quartic_real_center_radius':'min(a^(3/2),4/27-a), 0<a<4/27',
     'circle_polynomial':str(H),'boundary_squared_distance':str(distance),
     'global_geometric_argument':'see source: no cusp contact, no parallel smooth tangents, three double roots',
     'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
dest=Path(__file__).with_name('deltoid_inradius_verification.json')
dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'all_passed':True,'quartic_disk':out['quartic_unique_maximal_disk'],'output':str(dest)}))
