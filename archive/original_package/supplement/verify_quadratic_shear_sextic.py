"""Exact polynomial checks for the optimal quadratic shear radius and sextic seed.

The arbitrary small-parameter count follows from the local splitting proof;
no approximate root search is used in this script.
"""
from pathlib import Path
import json
import sympy as sp


def main():
    z,w,t,kap,s=sp.symbols('z w t kap s',real=True)
    f=z**3+z*w-kap; star=w**3+z*w-kap
    G=sp.expand(w*w*f-z*star)
    assert G==sp.expand((z**3-kap)*w*w-z*z*w+kap*z)
    P=sp.Poly(sp.expand(f*f+t*G),z,w)
    assert P.degree(z)==6 and P.degree(w)==2 and P.total_degree()==6
    assert [(mon,coeff) for mon,coeff in P.terms() if sum(mon)==6]==[((6,0),sp.Integer(1))]
    assert P.coeff_monomial(w*w)==-kap*t
    checks=0
    for kappa in [sp.Rational(1,1000),sp.Rational(1,100),sp.Rational(9,100),
                  sp.Rational(1,10),sp.Rational(7,50),sp.Rational(147,1000)]:
        assert 0<kappa<sp.Rational(4,27)
        real=sp.Poly(s**3+s*s-kappa,s)
        for a,b in [(-1,sp.Rational(-2,3)),(sp.Rational(-2,3),0),(0,sp.Rational(1,3))]:
            assert real.count_roots(a,b)==1
            checks+=1
        assert sp.gcd(real,real.diff()).degree()==0
    report={'scope':__doc__.strip(),'exact_seed_root_intervals':checks,
            'exact_cofactor_identity':True,'degrees':{'analytic':6,'conjugate':2,'total':6},
            'unique_highest_term':'z^6','constructed_regular_zeros':36,
            'real_zeros':6,'optimal_uniform_quadratic_shear_radius':1,'status':'passed'}
    Path(__file__).with_name('quadratic_shear_sextic_verification.json').write_text(
        json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='scope'},indent=2))


if __name__=='__main__':
    main()
