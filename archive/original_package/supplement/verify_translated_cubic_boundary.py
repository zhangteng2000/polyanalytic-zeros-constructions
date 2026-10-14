"""Exact audits for the all-order translated cubic-seed construction.
This script checks polynomial identities and degree support. The uniform
small-parameter and implicit-function arguments are supplied in the proof.
"""
from pathlib import Path
import json
import time
import sympy as s

U,T,A,H,C,D=s.symbols('U T A H C D')
z,w,kappa,tau=s.symbols('z w kappa tau')
f=z**3+z*w-kappa; g=w**3+z*w-kappa
AA=w*w; BB=-z; DD=-z*z*w; HH=-z**3; CC=-z*z*w
G=s.expand(AA*f+BB*g)
gb=s.groebner([f,g],z,w,order='grevlex',domain=s.QQ.frac_field(kappa))
def red(poly):return s.expand(gb.reduce(s.expand(poly))[1])

reports=[]
start=time.time()
for order in range(2,13):
    p=s.Rational(1,order)
    S=[s.Integer(3),U,U*U+2*p*T]
    for j in range(3,order+1):
        S.append(s.expand(U*S[-1]+p*T*S[-2]-p**3*C*S[-3]))
    Q=s.Poly(s.expand(S[order].subs({U:U+p*A,T:T-A*U+p*H},simultaneous=True)),U,T)
    assert Q.coeff_monomial(U**order)==1
    assert Q.coeff_monomial(U**(order-1))==0
    assert Q.coeff_monomial(U**(order-2)*T)==1
    global_poly=0
    term_reports=[]
    for (i,j),coef in Q.terms():
        weight=order-i-2*j
        assert weight>=0
        for (a,b,c),q in s.Poly(coef,A,H,C).terms():
            assert a+2*b+3*c==weight
        assert i+j<=order
        if j:
            # Replacing A by D changes every coefficient by a multiple of D-A.
            delta=s.expand(coef.subs(A,D)-coef)
            assert s.rem(delta,D-A,D)==0
            replacement=s.expand(coef.subs({A:DD,H:HH,C:CC}))
            term=s.Poly(s.expand(replacement*f**i*G**j),z,w)
            assert term.total_degree()<=3*order-j
            assert term.degree(w)<=order
        elif i==order:
            replacement=s.Integer(1)
            term=s.Poly(f**order,z,w)
        else:
            assert i<=order-2
            original=s.expand(coef.subs({A:AA,H:HH,C:CC}))
            replacement=red(original)
            assert red(replacement-original)==0
            term=s.Poly(s.expand(replacement*f**i),z,w)
            if replacement:
                assert term.total_degree()<=3*order-2
                assert term.degree(w)<=order
        global_poly+=tau**(order-i-j)*term.as_expr()
        term_reports.append([i,j,weight])
    pp=s.Poly(s.expand(global_poly),z,w)
    assert pp.total_degree()==3*order
    assert pp.degree(w)==order
    head=s.Add(*[q*z**i*w**j for (i,j),q in pp.terms() if i+j==3*order])
    assert s.expand(head-z**(3*order))==0
    assert s.expand(pp.as_expr().subs(tau,0)-f**order)==0
    row={'order':order,'local_terms':len(Q.terms()),'global_terms':len(pp.terms()),
         'total_degree':pp.total_degree(),'conjugate_degree':pp.degree(w),
         'highest_part':str(head),'all_checks_passed':True}
    reports.append(row)
    print(row,'elapsed',round(time.time()-start,2),flush=True)
out={'description':'Exact symbolic checks; no numerical zero-count search.',
     'orders':reports,'coefficient_weight_identity':'a+2*b+3*q=k-i-2*j',
     'normalized_parameter_defect':'abs((w^2+z^2*w)/z)',
     'all_passed':True}
Path(__file__).with_name('translated_cubic_boundary_verification.json').write_text(json.dumps(out,indent=2))
