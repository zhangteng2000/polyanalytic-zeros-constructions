"""Exact checks for the angular-separation theorem; Newton checks are diagnostic.

The proof of the infinite theorem and of the transcendental sector estimate
is in the accompanying TeX source, not inferred from these finite samples.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd, log, exp
import importlib.util
import json
import sympy as sp


def powers_of_endpoints(n, m, k):
    d=n-k
    return (F(d**d*(k+m)**(k-m),(n+m)**(n-m)),
            F(d**d*(k-m)**(k-m),(n-m)**(n-m)))


def rational_parameter(n, m, k):
    low,high=powers_of_endpoints(n,m,k)
    d=n-k
    mid=(exp(log(float(low))/d)+exp(log(float(high))/d))/2
    den=1
    while True:
        e=F(round(mid*den),den)
        if low<e**d<high:
            return e
        den*=2


def main():
    n,m,k,x,y=sp.symbols('n m k x y',real=True)
    d=n-k; K=k+m; L=n+m
    a=(k-m)/(n-m); b=(k+m)/(n+m)
    center=(a+b)/2; radius=(b-a)/2
    Dk=k*k+m*m-2*k*m*x
    Dn=n*n+m*m-2*n*m*x
    identities=[
        center**2-radius**2-a*b,
        (n*n-m*m)*(a+b)-2*(n*k-m*m),
        (n*n-m*m)*a*b-(k*k-m*m),
        k*k-2*k*n*x+n*n*(x*x+y*y)-m*m*(1-2*x+x*x+y*y)
        -(n*n-m*m)*((x-center)**2+y*y-radius**2),
        K*(m*m-k*m*x)*Dn-L*(m*m-n*m*x)*Dk
        -d*m*m*(1+x)*(-(k-m)*(n-m)+2*k*n*(1-x)),
        K*m*(k*x-m)-m*x*Dk-m*m*(1+x)*(2*k*x-K),
        L*(n-m*x)-Dn-m*(n-m)*(1+x),
        K*(k-m*x)-Dk-m*(k-m)*(1+x),
        4*k*k-(k+m)**2-(k-m)*(3*k+m),
        2*k*k-m*(k+m)-(k-m)*(2*k+m),
    ]
    assert all(sp.factor(q)==0 for q in identities)
    report={'scope':__doc__.strip(),'exact_symbolic_identities':len(identities)}

    checks=0; m2checks=0; broadchecks=0; witnesses=[]
    for nn in range(3,61):
        for mm in range(1,min(nn-1,11)):
            for kk in range(mm+1,nn):
                ell=gcd(nn-kk,kk+mm)
                sufficient=(mm<=2 or ell>=mm or 2*ell*(kk-mm)>=mm*mm)
                if not sufficient:
                    continue
                low,high=powers_of_endpoints(nn,mm,kk)
                assert low<high
                ep=rational_parameter(nn,mm,kk)
                assert low<ep**(nn-kk)<high
                assert ell==gcd(nn+mm,kk+mm)
                assert len({((kk+mm)*j)%(nn-kk) for j in range(nn-kk)})==(nn-kk)//ell
                assert nn+mm+1+2*ell<=nn*nn
                checks+=1
                m2checks+=(mm==2)
                broadchecks+=(mm>2)
                if nn<=10 and mm==2:
                    witnesses.append({'n':nn,'m':mm,'k':kk,'epsilon':str(ep),
                                      'maximum':nn+mm+1+2*ell})
    report['exact_interval_and_orbit_checks']=checks
    report['conjugate_degree_two_checks']=m2checks
    report['higher_conjugate_degree_checks']=broadchecks

    # The displayed degree-eight example is tested without floating point.
    lo,hi=powers_of_endpoints(8,2,3)
    assert lo==F(5**6,10**6) and hi==F(5**5,6**6)
    assert lo<F(4,7)**5<hi
    report['degree_eight_exact_parameter']={'epsilon':'4/7','count':21,
                                           'regular_nonzero_zeros':20,'origin_index':-2}

    spec=importlib.util.spec_from_file_location('v9',Path(__file__).with_name('verify_trinomial_global.py'))
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    samples=[(nn,2,kk) for nn in range(4,10) for kk in range(3,nn)]
    samples.extend([(9,3,8),(10,3,8),(12,3,8),(14,3,8),(16,4,6),(17,5,6)])
    diagnostics=[]
    for nn,mm,kk in samples:
        ep=rational_parameter(nn,mm,kk)
        ell=gcd(nn-kk,kk+mm)
        count,pos,neg=module.newton_count(nn,mm,kk,float(ep))
        predicted=nn+mm+1+2*ell
        assert count==predicted,(nn,mm,kk,ep,count,predicted)
        assert neg==ell and pos-neg==nn+mm,(nn,mm,kk,pos,neg)
        diagnostics.append({'n':nn,'m':mm,'k':kk,'epsilon':str(ep),
                            'count':count,'positive_nonzero':pos,'negative_nonzero':neg})
    report['rational_witnesses']=witnesses
    report['newton_diagnostics_not_certificates']=diagnostics
    report['status']='passed'
    Path(__file__).with_name('separated_harmonic_trinomial_verification.json').write_text(
        json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items()
                     if k not in ['scope','rational_witnesses','newton_diagnostics_not_certificates']},indent=2))
    print('Newton diagnostics:',len(diagnostics),'all passed')


if __name__=='__main__':
    main()
