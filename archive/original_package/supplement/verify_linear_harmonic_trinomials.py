"""Exact finite arithmetic checks and separate numerical diagnostics for v11.

The universal result follows from the logarithmic-convexity proof in the
TeX source. Newton searches below are diagnostics, not general certificates.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd, log, exp
import importlib.util
import json
import sympy as sp


def interval(n,k):
    d=n-k
    low=F(d**d*(k+1)**(k-1),(n+1)**(n-1))
    high=F(d**d*(k-1)**(k-1),(n-1)**(n-1))
    return low,high


def rational_maximizer(n,k):
    low,high=interval(n,k); d=n-k
    assert low<high
    mid=(exp(log(float(low))/d)+exp(log(float(high))/d))/2
    denominator=1
    while True:
        candidate=F(round(mid*denominator),denominator)
        if low<candidate**d<high:
            return candidate
        denominator*=2


def largest_proper_divisor(n):
    for p in range(2,int(n**.5)+1):
        if n%p==0:
            return n//p
    return 1


def main():
    n,k,x,y=sp.symbols('n k x y',real=True)
    d=n-k
    a=(k-1)/(n-1); b=(k+1)/(n+1)
    center=(a+b)/2; radius=(b-a)/2
    identities=[center**2-radius**2-a*b,
                (n*n-1)*(a+b)-2*(n*k-1),
                (n*n-1)*a*b-(k*k-1),
                (k*k-2*k*n*x+n*n*(x*x+y*y))
                -(1-2*x+x*x+y*y)
                -(n*n-1)*((x-center)**2+y*y-radius**2)]
    assert all(sp.factor(q)==0 for q in identities)
    report={'scope':__doc__.strip(),'symbolic_identities':len(identities)}

    checks=0; witnesses=[]
    for nn in range(3,76):
        for kk in range(2,nn):
            ep=rational_maximizer(nn,kk)
            low,high=interval(nn,kk)
            assert low<ep**(nn-kk)<high
            assert gcd(nn-kk,kk+1)==gcd(nn+1,kk+1)
            ell=gcd(nn-kk,kk+1)
            assert nn+2+2*ell<=3*nn-2
            checks+=1
            if nn<=10:
                witnesses.append({'n':nn,'k':kk,'epsilon':str(ep),
                                  'maximum':nn+2+2*ell})
    report['exact_rational_maximizer_checks']=checks

    optimization=0; equality=[]
    for nn in range(4,2001):
        values=[nn+2+2*gcd(nn+1,kk+1) for kk in range(2,nn)]
        endpoint=nn+4 if nn%2 else nn+2
        maximum=max(values+[nn+2,endpoint])
        asserted=nn+2+2*largest_proper_divisor(nn+1)
        assert maximum==asserted
        if maximum==3*nn-2:
            equality.append(nn)
        optimization+=len(values)
    assert equality==[5]
    report['exact_gcd_optimization_candidates']=optimization
    report['equality_degrees_from_4_to_2000']=equality

    # Independent radial elimination for the thirteen-zero quintic.
    r=sp.symbols('r')
    eps=sp.Rational(5,12)
    poly=sp.Poly(r**4-r+eps,r)
    assert poly.count_roots(0,sp.oo)==2
    assert eps**3>sp.Rational(1,16)
    assert eps**3<sp.Rational(27,256)
    report['exact_quintic_radial_certificate']={
        'epsilon':'5/12','positive_radial_roots':2,
        'ray_zeros':6,'off_ray_zeros':6,'origin':1,'total':13}

    spec=importlib.util.spec_from_file_location('v9',Path(__file__).with_name('verify_trinomial_global.py'))
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    diagnostics=[]
    for nn in range(3,11):
        for kk in range(2,nn):
            ep=rational_maximizer(nn,kk)
            predicted=nn+2+2*gcd(nn-kk,kk+1)
            count,positive,negative=module.newton_count(nn,1,kk,float(ep))
            assert count==predicted,(nn,kk,ep,count,predicted)
            assert negative==gcd(nn-kk,kk+1)
            # The common Newton helper counts orientations only away from zero.
            assert positive-negative==nn+1
            diagnostics.append({'n':nn,'k':kk,'epsilon':str(ep),
                                'predicted':predicted,'newton_count':count})
    report['rational_maximizer_examples']=witnesses
    report['newton_diagnostics_not_certificates']=diagnostics
    report['status']='passed'
    Path(__file__).with_name('linear_harmonic_trinomial_verification.json').write_text(
        json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:value for key,value in report.items()
                     if key not in ['scope','newton_diagnostics_not_certificates','rational_maximizer_examples']},indent=2))
    print('Newton diagnostics:',len(diagnostics),'all passed')


if __name__=='__main__':
    main()
