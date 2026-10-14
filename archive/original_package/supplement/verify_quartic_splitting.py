"""Exact symbolic diagnostics for the quartic splitting theorem.

Requires Python 3 and SymPy. These checks supplement the all-r proof;
they do not claim a certified numerical choice of the small parameters.
"""
from pathlib import Path
import json
import sympy as s

if not __debug__:
    raise RuntimeError('Run without Python optimization: assertions are required.')

u, v, a, ac, b, bc = s.symbols('u v a ac b bc')
H = lambda x, y: x**4+b*x*x*y+b*b*y*y/8-b*b*bc*x/16+bc**4/256
G = a*u+b*v
T1 = u*u*G
T2 = G**2/8+a*u*G/4+b*ac*u*u/4
T3 = a*a*G/16+b*ac*(a*u+G)/16-b*b*bc*u/16
T4 = a**4/256+b*a*a*ac/64+b*b*ac*ac/128-b*b*bc*a/64+bc**4/256
assert s.expand(H(u+a/4,v+ac/4)-(u**4+T1+T2+T3+T4)) == 0
print('Local translation identity passed.',flush=True)

z,w = s.symbols('z w')
lam=1+2*s.I
domain=s.QQ_I
def poly(expr):
    return s.Poly(expr,z,w,domain=domain)

records=[]
for r in range(3,25):
    print('Checking r =',r,flush=True)
    N=2*r+1
    # epsilon=1 is used only for a degree identity, never as a zero-count example.
    Q=poly(1)
    for j in range(1,r+1): Q*=poly(z*w-j*j)
    c1=-sum(j*j for j in range(1,r+1))
    c2=sum(j*j*k*k for j in range(1,r+1) for k in range(j+1,r+1))
    assert Q.coeff_monomial(z**(r-1)*w**(r-1))==c1
    assert Q.coeff_monomial(z**(r-2)*w**(r-2))==c2
    f=poly(z**N)+Q; fs=poly(w**N)+Q
    A=poly(lam*w**(r+1)); B=poly(-lam*z**r)
    As=poly(s.conjugate(lam)*z**(r+1)); Bs=poly(-s.conjugate(lam)*w**r)
    g=A*f+B*fs
    d2=poly(lam**2*c1*z**(r-1)*w**r/4)
    d3=poly(lam**3*(c1*z**(r-1)*w**(2*r+1)+c2*z**(r-2)*w**(2*r))/16)
    terms=[
        f*f*g,
        g*g/s.Integer(8)+A*f*g/s.Integer(4)+B*As*f*f/s.Integer(4)-d2*f*fs,
        A*A*g/s.Integer(16)+B*As*(A*f+g)/s.Integer(16)-B*B*Bs*f/s.Integer(16)-d3*fs,
        A**4/s.Integer(256)+B*A*A*As/s.Integer(64)+B*B*As*As/s.Integer(128)
        -B*B*Bs*A/s.Integer(64)+Bs**4/s.Integer(256)-poly(lam**4*w**(2*r+3)/256)*fs,
    ]
    terms=[poly(p.as_expr()) for p in terms]
    bounds=[7*r+4,6*r+4,5*r+4,4*r+4]
    assert g.degree(w)<=2*r
    assert g.coeff_monomial(w**(2*r)*z**(r-1))==lam*c1
    assert g.coeff_monomial(w**(2*r-1)*z**(r-2))==lam*c2
    for p,d in zip(terms,bounds):
        assert p.degree(w)<=4*r
        assert p.total_degree()<=d<4*N
    assert (f**4).coeff_monomial(z**(4*r)*w**(4*r))==1
    assert (f**4).coeff_monomial(z**(4*N))==1
    records.append({'r':r,'m':4*r,'n':4*N,
                    'correction_w_degrees':[p.degree(w) for p in terms],
                    'correction_total_degrees':[p.total_degree() for p in terms],
                    'E':(4*N)**2,'M':(4*r)**2+(4*N)**2})

# A generic radial polynomial tests cancellations without numeric coefficients.
r=3; N=7
q0,q1,q2=s.symbols('q0 q1 q2')
Q=(z*w)**3+q2*(z*w)**2+q1*z*w+q0
f=z**N+Q; fs=w**N+Q
A=lam*w**4; B=-lam*z**3
As=s.conjugate(lam)*z**4; Bs=-s.conjugate(lam)*w**3
g=s.expand(A*f+B*fs)
d2=lam**2*q2*z**2*w**3/4
d3=lam**3*(q2*z**2*w**7+q1*z*w**6)/16
terms=[f*f*g,g*g/8+A*f*g/4+B*As*f*f/4-d2*f*fs,
       A*A*g/16+B*As*(A*f+g)/16-B*B*Bs*f/16-d3*fs,
       A**4/256+B*A*A*As/64+B*B*As*As/128-B*B*Bs*A/64+Bs**4/256-lam**4*w**9*fs/256]
for p in terms:
    pp=s.Poly(s.expand(p),z,w)
    assert pp.degree(w)<=12 and pp.total_degree()<28

rho=s.simplify(lam/s.conjugate(lam))
assert s.simplify(rho+1/rho)==-s.Rational(6,5)
eps,S,T=s.symbols('epsilon S T',positive=True)
assert s.simplify(eps**2*S/(4*(3*eps/4))-eps*S/3)==0
assert s.simplify((eps**2*S+eps**4*T/(3*eps/4)**2)/16-eps**2*(S/16+T/9))==0

report={'status':'passed','arithmetic':'exact Gaussian rational and symbolic arithmetic',
        'local_translation_identity':True,'generic_radial_identity_r3':True,
        'degree_cases':len(records),'degree_range':'3 <= r <= 24',
        'degree_specialization':'epsilon=1, for support identities only',
        'normalized_error_bound_identities':2,
        'phase_ratio_trace':'-6/5, a noninteger rational',
        'scope':'The infinite sharpness theorem is proved analytically. No fixed numerical small parameter or floating-point zero count is certified here.',
        'cases':records}
target=Path(__file__).with_name('quartic_splitting_verification.json')
target.write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='cases'},indent=2))
