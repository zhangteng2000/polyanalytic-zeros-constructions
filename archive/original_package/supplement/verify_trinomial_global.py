"""Exact identities and separately labelled numerical diagnostics for v9.

The universal zero counts follow from the analytic proof in the TeX source.
Newton samples below are diagnostics, not certificates for all degrees.
"""
from pathlib import Path
from math import gcd, pi, cos, atan2, ceil
import json
import numpy as np
import sympy as sp


def phase(t, n, m, k):
    d = n-k
    arg_num = t+pi+atan2((k/m)*np.sin(t), 1-(k/m)*cos(t))
    arg_den = atan2(-m*np.sin(t), n-m*cos(t))
    return ((m+k)*arg_num-(n+m)*arg_den)/d


def critical_radii(alpha, n, m, k):
    d, K = n-k, m+k
    start = K*pi/d
    first = ceil((start-alpha)/(2*pi))
    radii = []
    for j in range(first, first+K):
        target = alpha+2*pi*j
        lo, hi = 0., 2*pi*d
        for _ in range(65):
            mid = (lo+hi)/2
            if phase(mid,n,m,k)<target:
                lo=mid
            else:
                hi=mid
        t=(lo+hi)/2
        radii.append(d/(abs(n-m*np.exp(1j*t))**((n-m)/d)
                         *abs(k-m*np.exp(1j*t))**((m-k)/d)))
    return np.array(radii)


def newton_count(n,m,k,epsilon):
    radii=np.array([.15,.3,.5,.75,1.,1.4,2.,3.,4.])
    angles=np.arange(480)*2*pi/480
    z=(radii[:,None]*np.exp(1j*angles)).ravel()
    for _ in range(95):
        F=z**n-z**k+epsilon*z.conj()**m
        A=n*z**(n-1)-k*z**(k-1)
        B=m*epsilon*z.conj()**(m-1)
        J=abs(A)**2-abs(B)**2
        J=np.where(abs(J)>1e-80,J,1e-80)
        step=(A.conj()*F-B*F.conj())/J
        step/=np.maximum(1,abs(step)/2)
        z-=step
    F=z**n-z**k+epsilon*z.conj()**m
    scale=abs(z)**n+abs(z)**k+abs(epsilon)*abs(z)**m
    # Absolute residuals alone would falsely accept small points near
    # the deliberately multiple origin when k is large.
    good=(abs(F)<1e-10*np.maximum(scale,1e-200))&(abs(z)>.05)&(abs(z)<5)
    rounded=np.unique(np.round(z[good].real,7)+1j*np.round(z[good].imag,7))
    roots=[]
    for a in rounded:
        if all(abs(a-b)>2e-6 for b in roots):
            roots.append(a)
    roots=np.array(roots)
    A=n*roots**(n-1)-k*roots**(k-1)
    B=m*epsilon*roots.conj()**(m-1)
    J=abs(A)**2-abs(B)**2
    assert np.all(abs(J)>1e-7*(abs(A)**2+abs(B)**2))
    return len(roots)+1, int(sum(J>0)), int(sum(J<0))


def main():
    n,m,k,x,v=sp.symbols('n m k x v', real=True)
    d=n-k; K=m+k; L=n+m
    numerator=sp.expand(
        K*(m*m-k*m*x)*(n*n+m*m-2*n*m*x)
        -L*(m*m-n*m*x)*(k*k+m*m-2*k*m*x))
    factored=d*m*m*(1+x)*((m-k)*(n-m)+2*k*n*(1-x))
    identities=[
        numerator-factored,
        m*m*d*d-(n*k-m*m)**2-(m*m-k*k)*(n*n-m*m),
        (d+1+2*K)-(n+2*m+k+1),
        (L-K)-d,
    ]
    assert all(sp.expand(q)==0 for q in identities)
    report={'scope':__doc__.strip(),'symbolic_identities':len(identities)}

    # Exact rational positivity in the interior of every tested cosine
    # interval; the general factorization proves it beyond this range.
    positivity=0
    for nn in range(3,31):
        for mm in range(2,nn):
            for kk in range(1,mm):
                for xx in [sp.Rational(-1),sp.Rational(-3,4),sp.Rational(0),
                           sp.Rational(3,4),sp.Rational(1)]:
                    A=(mm-kk)*(nn-mm)+2*kk*nn*(1-xx)
                    assert A>0
                    assert kk*kk+mm*mm-2*kk*mm*xx>0
                    assert nn*nn+mm*mm-2*nn*mm*xx>0
                    positivity+=1
    report['exact_argument_positivity_checks']=positivity

    orbit_checks=0
    for nn in range(3,101):
        for mm in range(2,nn):
            for kk in range(1,mm):
                dd,KK=nn-kk,mm+kk
                ell=gcd(dd,KK)
                assert len({(KK*j)%dd for j in range(dd)})==dd//ell
                assert (2*KK)%ell==0
                assert ((nn+2*mm+kk+1)-(dd+1))//ell==2*KK//ell
                orbit_checks+=1
    report['exact_symmetry_orbit_checks']=orbit_checks

    # An independent elimination and rational Sturm isolation for
    # F=z^3-z+3*bar(z)^2: three y=0 zeros and six y!=0 zeros.
    X,Y=sp.symbols('X Y',real=True)
    z=X+sp.I*Y
    f=sp.expand(z**3-z+3*sp.conjugate(z)**2)
    real,imag=sp.re(f),sp.im(f)
    assert sp.expand(imag-Y*(3*X*X-Y*Y-1-6*X))==0
    y2=3*X*X-6*X-1
    poly=sp.Poly(sp.expand(real).subs(Y**2,y2).expand(),X)
    assert poly.as_expr()==-8*X**3+12*X**2+20*X+3
    intervals=sp.polys.polytools.intervals(poly,eps=sp.Rational(1,10**8))
    assert len(intervals)==3
    for ((lo,hi),mult) in intervals:
        assert mult==1 and not (lo<=1<=hi)
        assert min(y2.subs(X,lo),y2.subs(X,hi))>0
    assert sp.Poly(real.subs(Y,0),X).count_roots(-sp.oo,sp.oo)==3
    report['exact_cubic_elimination']={'nonzero_y_roots':6,'real_axis_roots':3,'total':9}

    samples=[]
    for nn,mm,kk in [(3,2,1),(4,3,2),(5,4,3),(7,5,3),(9,7,5),(12,11,10)]:
        numerical,positive,negative=newton_count(nn,mm,kk,3.)
        expected=nn+2*mm+kk+1
        assert numerical==expected
        assert (positive,negative)==(nn+mm,mm+kk)
        samples.append({'degrees':[nn,mm,kk],'epsilon':[3,0],
                        'predicted':expected,'newton_count':numerical})
    for nn,mm,kk in [(5,3,1),(7,5,3),(8,5,2)]:
        alpha=pi/(2*(nn-kk))
        rs=critical_radii(alpha,nn,mm,kk)
        unique_radii=[]
        for r in sorted(rs):
            if not unique_radii or abs(r-unique_radii[-1])>1e-8:
                unique_radii.append(r)
        ell=gcd(nn-kk,mm+kk)
        assert len(unique_radii)==(mm+kk)//ell
        edges=[0.]+unique_radii+[3.]
        for lo,hi in zip(edges[:-1],edges[1:]):
            r=(lo+hi)/2
            epsilon=r*np.exp(1j*alpha)
            q=int(sum(rs<r))
            expected=nn-kk+1+2*q
            numerical,positive,negative=newton_count(nn,mm,kk,epsilon)
            assert numerical==expected,(nn,mm,kk,epsilon,numerical,expected)
            assert (positive,negative)==(nn-kk+q,q)
            samples.append({'degrees':[nn,mm,kk],
                            'epsilon':[float(epsilon.real),float(epsilon.imag)],
                            'predicted':expected,'newton_count':numerical})
    report['numerical_diagnostics_not_certificates']=samples
    report['numerical_sample_count']=len(samples)
    report['status']='passed'
    out=Path(__file__).with_name('trinomial_global_verification.json')
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='numerical_diagnostics_not_certificates'},indent=2))


if __name__=='__main__':main()
