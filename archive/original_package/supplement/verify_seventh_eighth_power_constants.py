"""Exact rational certificates for the optimal unit-triple constants C7,C8.
The accompanying proof reduces the global problem to these identities and
positive Bernstein coefficients. No floating comparison is accepted.
"""
from pathlib import Path
import hashlib,json
import sympy as s
if not __debug__:
    raise RuntimeError('Run without -O: proof assertions must remain enabled')
U,V,R,t,v,x=s.symbols('U V R t v x')
P=[s.S(3),U,U**2-2*V]
for k in range(3,9):P.append(s.expand(U*P[-1]-V*P[-2]+P[-3]))
A7=7*(1-R**2)*(1-2*R**2)+7*R**3*(1-R**2)*v+R**6*v**2+7*R*(1-R**2)/v
A8=-8*(1-R**2)*(1-2*R**2)+4*R*(1-R**2)*(3-5*R**2)*v+8*R**4*(1-R**2)*v**2+R**7*v**3+2*R**3/v
for k,A,shift in [(7,A7,0),(8,A8,2)]:
    assert s.expand(P[k].subs({U:R*t,V:R/t})/(R*t)*t**shift-A.subs(v,t**3))==0
f7=R**6-7*R**5+14*R**4-21*R**2+7*R+7
f8=R**7+8*R**6+20*R**5+8*R**4-30*R**3-24*R**2+12*R+8
assert s.expand(A7.subs(v,1)-f7)==0
assert s.expand(-A8.subs(v,-1)-f8)==0
L=s.expand(A8*A8.subs(v,1/v))
H=L.coeff(v,0)
for j in range(1,5):
    assert s.expand(L.coeff(v,j)-L.coeff(v,-j))==0
    H+=2*L.coeff(v,j)*s.chebyshevt(j,x)
H=s.expand(H)
assert s.expand(H.subs(x,-1)-f8**2)==0
D=s.cancel(-s.diff(H,x)/R)
assert s.Poly(D,R,x).degree_list()==(12,3)
def bernstein(poly,variables,intervals,degrees):
    q=s.Poly(s.expand(poly.subs({z:a+(b-a)*z for z,(a,b) in zip(variables,intervals)},simultaneous=True)),*variables)
    assert all(q.degree(z)<=n for z,n in zip(variables,degrees))
    import itertools
    coeff={}
    for idx in itertools.product(*(range(n+1) for n in degrees)):
        b=0
        for power in itertools.product(*(range(j+1) for j in idx)):
            term=q.coeff_monomial(s.prod(z**h for z,h in zip(variables,power)))
            term*=s.prod(s.binomial(j,h)/s.binomial(n,h) for j,n,h in zip(idx,degrees,power))
            b+=term
        coeff[idx]=s.factor(b)
    reconstructed=sum(b*s.prod(s.binomial(n,j)*z**j*(1-z)**(n-j) for z,n,j in zip(variables,degrees,idx)) for idx,b in coeff.items())
    assert s.expand(reconstructed-q.as_expr())==0
    return coeff
bc=bernstein(D,[R,x],[(s.S(0),s.Rational(3,8)),(-s.S(1),s.S(1))],[12,3])
assert min(bc.values())==s.Rational(303335857569,4294967296)>0
rows=[]
for k,f,right,expected in [(7,f7,s.Rational(3,7),s.Rational(50784,2401)),(8,f8,s.Rational(3,8),s.S(48))]:
    concavity=bernstein(-s.diff(f,R,2),[R],[(s.S(0),right)],[k-3])
    assert min(concavity.values())==expected>0
    derivative=s.diff(f,R)
    assert derivative.subs(R,0)>0>derivative.subs(R,right)
    assert f.subs(R,0)==k
    lo,hi=s.S(0),right
    for step in range(120):
        mid=(lo+hi)/2
        if derivative.subs(R,mid)>0:lo=mid
        else:hi=mid
    assert derivative.subs(R,lo)>0>derivative.subs(R,hi)
    # Rational interval arithmetic evaluates the polynomial at the unique
    # critical point. Negative monomials reverse the endpoint order.
    flo=fhi=s.S(0)
    for (power,),coefficient in s.Poly(f,R).terms():
        a,b=lo**power,hi**power
        if coefficient<0:a,b=b,a
        flo+=coefficient*a;fhi+=coefficient*b
    assert flo>k
    den=10**12;clo=s.floor(flo*den)/den;chi=s.ceiling(fhi*den)/den
    assert clo<=flo<fhi<=chi and chi-clo==s.Rational(1,den)
    rows.append({'order':k,'radial_polynomial':str(f),'stationary_polynomial':str(derivative),
                 'root_interval':[str(lo),str(hi)],'constant_interval':[str(clo),str(chi)],
                 'constant_decimal_for_reading':[str(s.N(clo,16)),str(s.N(chi,16))],
                 'concavity_bernstein_coefficients':[str(a) for a in concavity.values()],
                 'concavity_minimum':str(expected),
                 'equality_traces':('rho * omega^j' if k==7 else '-rho * omega^j')})
assert 1-2*s.Rational(3,7)**2>0
assert 8-24*s.Rational(3,8)**2-30*s.Rational(3,8)**3>0
out={'all_passed':True,'arithmetic':'exact rational polynomial arithmetic',
     'orders':rows,'angular_polynomial':str(D),
     'angular_bernstein_degree':[12,3],
     'angular_bernstein_coefficients':[[str(bc[i,j]) for j in range(4)] for i in range(13)],
     'angular_bernstein_minimum':str(min(bc.values())),
     'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
dest=Path(__file__).with_name('seventh_eighth_power_constants_verification.json')
dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'all_passed':True,'intervals':[r['constant_decimal_for_reading'] for r in rows],'output':str(dest)}))
