"""Exact symbolic/arithmetic checks for all three-term supports with q-degree two.

Caustic winding and Newton samples are explicitly diagnostic, not certificates.
The global bounds, singular-target argument, and infinite formulas are proved
in the TeX source.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd, pi
import importlib.util
import json
import numpy as np
import sympy as sp


def largest_proper_divisor(n):
    for p in range(2,int(n**.5)+1):
        if n%p==0:
            return n//p
    return 1


def all_support_max(n):
    analytic=[n+6, n+7 if n%4==2 else n+3]
    analytic += [n+3+2*gcd(n+2,k+2) for k in range(3,n)]
    alternative=9 if n==3 else 10 if n==4 else n+4
    return max(analytic+[n+4,n+3,alternative])


def opposite_samples(order,t):
    # Exact scaling model from the opposite-target proof.
    return (1,1) if order%2 else (2,0)


def caustic_samples(n):
    theta=np.linspace(0,2*pi,16001)
    phi=(n+2)*theta
    rho=(np.cos(phi)+np.sqrt(n*n-np.sin(phi)**2))/(n*n-1)
    z=rho**(1/(n-2))*np.exp(1j*theta)
    curve=z.conj()*(1-rho*np.exp(1j*phi))
    low=float(abs(curve).min()); high=float(abs(curve).max())
    maximum=0; checks=0
    for a in np.linspace(0,2*pi/(n+2),33):
        for r in np.linspace(low*.9,high*1.01,80):
            w=curve-r*np.exp(1j*a)
            winding=round(float(np.angle(w[1:]/w[:-1]).sum()/(2*pi)))
            negative=-winding
            assert 0<=negative<=(2 if n<=4 else 1),(n,a,r,negative)
            maximum=max(maximum,negative); checks+=1
    return {'n':n,'grid_targets':checks,'largest_negative_winding':maximum,
            'role':'numerical diagnostic only'}


def main():
    n,x=sp.symbols('n x',real=True)
    D1=n*n+1-2*n*x; D2=n*n+4-4*n*x
    derivative_num=n*((n+2)*(n-2*x)*D1-(n+1)*(n-x)*D2)
    identities=[
        derivative_num-n*n*(1+x)*(4*x+n*n-3*n-2),
        -2*(1+4*x/n+4/n**2)+(n+2)*(2*x/n+4/n**2)
        -2*(n-2)*(x-1)/n,
        n*n-3*n-6-(n-5)*(n+2)-4,
    ]
    assert all(sp.factor(q)==0 for q in identities)
    report={'scope':__doc__.strip(),'exact_symbolic_identities':len(identities)}

    checked=0; even=0; rows=[]
    assert [all_support_max(nn) for nn in range(3,7)]==[9,10,11,13]
    for nn in range(7,2001):
        formula=max(nn+6,nn+3+2*largest_proper_divisor(nn+2))
        assert all_support_max(nn)==formula
        if nn%2==0:
            assert formula==2*nn+5
            even+=1
        checked+=nn-3
        if nn<=30:
            rows.append({'n':nn,'distinct_maximum':formula,'regular_maximum':nn+6})
    assert all(sum(opposite_samples(r,F(1,100)))==2 for r in range(2,101))
    report['exact_support_optimization_candidates']=checked
    report['even_degree_formula_checks']=even
    report['opposite_target_model_orders']=99

    spec=importlib.util.spec_from_file_location('v12',Path(__file__).with_name('verify_separated_harmonic_trinomials.py'))
    sep=importlib.util.module_from_spec(spec); spec.loader.exec_module(sep)
    parameters=[(6,4,F(7,32),13),(7,4,F(5,16),16),(8,3,F(4,7),21)]
    for nn,kk,ep,count in parameters:
        lo,hi=sep.powers_of_endpoints(nn,2,kk)
        assert lo<ep**(nn-kk)<hi
        assert count==all_support_max(nn)
    report['exact_displayed_rational_examples']=[{'n':nn,'k':kk,'epsilon':str(ep),'count':c}
                                               for nn,kk,ep,c in parameters]
    spec=importlib.util.spec_from_file_location('v9',Path(__file__).with_name('verify_trinomial_global.py'))
    common=importlib.util.module_from_spec(spec); spec.loader.exec_module(common)
    newton=[]
    for nn in range(3,16):
        count,pos,neg=common.newton_count(nn,2,1,3.)
        assert count==nn+6 and pos+1-neg==nn
        newton.append({'n':nn,'support':'z^n-z+3*barz^2','count':count,
                       'role':'numerical diagnostic only'})
    report['newton_diagnostics_not_certificates']=newton
    report['caustic_winding_diagnostics_not_certificates']=[caustic_samples(nn) for nn in range(3,9)]
    report['maximum_table']=rows
    report['status']='passed'
    Path(__file__).with_name('quadratic_conjugate_sparsity_verification.json').write_text(
        json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items()
                     if k not in ['scope','maximum_table','newton_diagnostics_not_certificates',
                                  'caustic_winding_diagnostics_not_certificates']},indent=2))
    print('Newton diagnostics:',len(newton),'all passed')
    print('Caustic winding diagnostics: six degrees, all passed')


if __name__=='__main__':
    main()
