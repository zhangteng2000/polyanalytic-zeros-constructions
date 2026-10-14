"""Exact finite checks for the k=m endpoint, with labelled Newton diagnostics.

Universal assertions are proved in polyanalytic_sharp_v3.tex, not inferred
from this finite enumeration or the numerical samples.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd, floor, ceil, pi, atan2, cos, sin
import importlib.util
import json
import numpy as np
import sympy as sp


def unit_count(n, m, alpha_over_pi):
    """Count positive signs exactly; no trigonometric floating point."""
    d=n-m
    offset=(m-d*alpha_over_pi)/(2*n)
    return 1+sum(0<(offset+F(n+m,n)*j)%2<1 for j in range(2*n))


def arc_coverage(n, m, alpha_over_2pi):
    """Count lifted open argument intervals from all d domain arcs."""
    d=n-m
    half=F(m,2*d)
    count=0
    for j in range(d):
        center=F(2*m*j,d)
        lo=center-half-alpha_over_2pi
        hi=center+half-alpha_over_2pi
        count+=ceil(hi)-floor(lo)-1
    return count


def critical_radii(n,m,alpha):
    d=n-m
    def phase(t,j):
        return (m*(t-pi)-(n+m)*atan2(-m*sin(t),n-m*cos(t))
                +4*pi*m*j)/d
    roots=[]
    for j in range(d):
        lower=phase(0.,j); upper=phase(2*pi,j)
        for a in range(floor((lower-alpha)/(2*pi))+1,
                       ceil((upper-alpha)/(2*pi))):
            target=alpha+2*pi*a
            lo,hi=0.,2*pi
            for _ in range(62):
                mid=(lo+hi)/2
                if phase(mid,j)<target:
                    lo=mid
                else:
                    hi=mid
            t=(lo+hi)/2
            roots.append(d/abs(n-m*np.exp(1j*t)))
    return np.array(roots)


def main():
    n,m,x=sp.symbols('n m x',real=True)
    den=n*n+m*m-2*n*m*x
    identities=[m*den-(n+m)*(m*m-m*n*x)-(n-m)*m*n*(1+x),
                (n*n-m*m)-m*(n-m)-n*(n-m),
                (n*n)-(n*m)-n*(n-m)]
    assert all(sp.expand(z)==0 for z in identities)
    report={'scope':__doc__.strip(),'symbolic_identities':len(identities)}

    unit_tests=0; orbit_tests=0; coverage_tests=0; real_tests=0
    for nn in range(2,76):
        for mm in range(1,nn):
            dd=nn-mm; ell=gcd(dd,2*mm); J=2*mm//ell; D=dd//ell
            assert gcd(nn+mm,2*nn)==ell
            assert (2*nn//ell)%2==J%2
            assert J%2==0 or ell%2==0
            orbit_tests+=1
            allowed=({nn-ell+1,nn+1} if J%2==0
                     else {nn+1-ell//2,nn+1+ell//2})
            g=gcd(dd,mm)
            special=F(mm,dd)
            # Includes endpoint alignment, either side, and independent phases.
            angles=[F(0),F(1),F(1,7),F(13,23),special,
                    special+F(g,2*dd),special-F(g,2*dd),
                    special+F(g,dd)]
            seen={unit_count(nn,mm,A) for A in angles}
            assert seen==allowed,(nn,mm,seen,allowed)
            assert unit_count(nn,mm,special)==min(allowed)
            unit_tests+=len(angles)
            Imax=ell*((J+1)//2)
            generic=F(1,4*dd)
            other=generic+F(1,2*D)
            counts=[arc_coverage(nn,mm,A) for A in [generic,other]]
            expected={ell*(J//2),ell*((J+1)//2)}
            assert set(counts)==expected,(nn,mm,counts,expected)
            for A in [F(0),F(1,2),F(mm,2*dd),generic,other]:
                C=arc_coverage(nn,mm,A)
                assert 0<=C<=Imax
                coverage_tests+=1
            if J%2:
                Ireal=max(arc_coverage(nn,mm,F(0)),
                          arc_coverage(nn,mm,F(1,2)))
                exceptional=J%4==3 and D%2==0
                computed=max(dd+1+2*Ireal,nn+mm+1,max(allowed))
                predicted=nn+mm+1 if exceptional else nn+mm+ell+1
                assert computed==predicted
                real_tests+=1
    report.update(exact_unit_circle_counts=unit_tests,
                  exact_orbit_and_parity_checks=orbit_tests,
                  exact_arc_coverage_checks=coverage_tests,
                  exact_real_parameter_maximum_checks=real_tests)

    # Verify the explicit unit-circle roots separately by substitution.
    substitutions=0
    for nn in range(2,19):
        for mm in range(1,nn):
            for A in [F(0),F(1),F(1,7),F(mm,nn-mm)]:
                alpha=float(A)*pi
                theta=(alpha+(2*np.arange(2*nn)+1)*pi)/(2*nn)
                c=2*(-1.)**np.arange(2*nn)*np.sin(mm*theta-alpha/2)
                roots=c[c>1e-10]**(1/(nn-mm))*np.exp(1j*theta[c>1e-10])
                value=roots**nn-roots**mm+np.exp(1j*alpha)*roots.conj()**mm
                scale=abs(roots)**mm*(2+abs(roots)**(nn-mm))
                assert np.all(abs(value)<2e-12*scale)
                assert len(roots)+1==unit_count(nn,mm,A)
                substitutions+=1
    report['numerical_unit_root_substitutions']=substitutions

    spec=importlib.util.spec_from_file_location('v9',Path(__file__).with_name('verify_trinomial_global.py'))
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    diagnostics=[]
    cases=[(2,1,.6),(3,1,.75),(4,1,.8),(5,3,-.99),
           (7,3,.99j),(7,3,.99),(7,3,-.99),
           (7,3,2.),(7,5,.99),(7,5,-.99),(9,5,.99j)]
    for nn,mm,ep in cases:
        if abs(ep)>1:
            predicted=nn+mm+1
        else:
            radii=critical_radii(nn,mm,np.angle(ep))
            assert np.all(abs(radii-abs(ep))>1e-5)
            predicted=nn-mm+1+2*int(sum(radii<abs(ep)))
        found,positive,negative=module.newton_count(nn,mm,mm,ep)
        assert found==predicted,(nn,mm,ep,found,predicted)
        diagnostics.append({'n':nn,'m':mm,'epsilon':[float(np.real(ep)),float(np.imag(ep))],
                            'predicted_distinct':predicted,'newton_distinct':found,
                            'positive_nonzero':positive,'negative_nonzero':negative})
    report['newton_diagnostics_not_certificates']=diagnostics
    report['status']='passed'
    target=Path(__file__).with_name('trinomial_endpoint_verification.json')
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['scope','newton_diagnostics_not_certificates']},indent=2))
    print('Newton diagnostics:',len(diagnostics),'all passed')


if __name__=='__main__':
    main()
