"""Exact algebraic checks for pure-head contact splitting.

The arbitrary-degree sharp endpoint additionally uses Geyer's critically
fixed construction, as cited and proved by inversion in the TeX source.
"""
from pathlib import Path
import json
import sympy as sp


def main():
    n,m,q=sp.symbols('n m q',integer=True)
    assert sp.expand((4*m+1)*(n+1)-2*m*(m+2)-(2*q+1)
                     -((4*m+1)*n-2*m*m-2*q))==0
    assert sp.expand((3*n-2)+(2*n+1)-(5*(n+1)-6))==0
    checks=0
    for nn in range(2,151):
        for mm in range(1,nn):
            for qq in [2,nn]:
                general=(4*mm+1)*nn-2*mm*(mm+2)
                contact=(4*mm+1)*nn-2*mm*mm-2*qq
                assert general-contact==2*(qq-2*mm)
                assert (nn+1)**2-(2*qq+1)>=nn**2
                checks+=1
        assert min(nn*nn,3*nn-2)==3*nn-2
    assert min(6**2,9*6-2*2**2-2*6)==34
    assert min(10**2,13*10-2*3**2-2*10)==92

    r=sp.symbols('r')
    eps=sp.Rational(1,10)
    positive=sp.Poly(r**3-r**2+eps,r)
    negative=sp.Poly(r**3+r**2-eps,r)
    assert eps<sp.Rational(4,27)
    assert positive.count_roots(0,sp.oo)==2
    assert negative.count_roots(0,sp.oo)==1
    assert sp.gcd(positive,positive.diff()).degree()==0
    assert sp.gcd(negative,negative.diff()).degree()==0
    assert positive.eval(sp.Rational(2,3))!=0

    u,v=sp.symbols('u v',nonzero=True)
    inversion=0
    for nn in range(2,25):
        coeff=sp.symbols('c2:'+str(nn+1))
        g=sum(coeff[j-2]*u**j for j in range(2,nn+1))
        A=sp.expand(u**nn*g.subs(u,1/u))
        assert sp.Poly(A,u).degree()==nn-2
        assert A.subs(u,0)==coeff[-1]
        assert sp.expand(v*A-u**nn-u**nn*v*(g.subs(u,1/u)-1/v))==0
        inversion+=1
    report={'scope':__doc__.strip(),'exact_bound_comparisons':checks,
            'exact_inversion_identities':inversion,
            'sixth_order_pure_contact_bound_at_2_6':34,
            'tenth_order_pure_contact_bound_at_3_10':92,
            'exact_cubic_example':{'polynomial':'z*barz-1/10-z^3',
                'positive_phase_radial_roots':2,'negative_phase_radial_roots':1,
                'roots_per_radius':3,'total_regular_zeros':9},'status':'passed'}
    Path(__file__).with_name('pure_contact_verification.json').write_text(
        json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='scope'},indent=2))


if __name__=='__main__':
    main()
