"""Finite exact checks for local fiber counts and real-coefficient extremizers.

Requires Python 3 and SymPy. Run without -O. This script checks algebra and
finite orbit counts; it does not replace the all-degree proofs in the TeX.
"""
from pathlib import Path
from fractions import Fraction
from itertools import combinations_with_replacement
from math import gcd
import json
import time
import sympy as S

if not __debug__:
    raise RuntimeError('Run this verifier without Python optimization (-O).')

u,v,z,w,t = S.symbols('u v z w t')


def poly_T(k):
    p = lambda x: S.Poly(x,u,v,domain=S.QQ)
    rec = [p(3),p(u),p(u*u-2*v)]
    for j in range(3,k+1):
        rec.append(p(u)*rec[-1]-p(v)*rec[-2]+rec[-3])
    return rec[k]


def swap(P):
    return S.Poly.from_dict({(b,a):c for (a,b),c in P.terms()},(u,v),domain=S.QQ)


def discriminant_checks(k):
    p = lambda x: S.Poly(x,u,v,domain=S.QQ)
    T = poly_T(k); R = swap(T)
    J = T.diff(u)*R.diff(v)-T.diff(v)*R.diff(u)
    D = p(u*u*v*v-4*u**3-4*v**3+18*u*v-27)
    Dout = T*T*R*R-4*T**3-4*R**3+18*T*R-p(27)
    assert (k**4*Dout-D*J*J).is_zero
    assert all(a+b<k for (a,b),c in (T-p(u**k)).terms())
    return {'k':k, 'terms':len(T.terms()), 'jacobian_identity':True}


def cusp_orbits(k,j):
    # lambda^k=omega^j, product of the three lambda's is one.
    # Exponents are exact rational turns, reduced modulo one.
    multisets=set()
    for a in range(k):
        x=Fraction(j+3*a,3*k)%1
        for b in range(k):
            y=Fraction(j+3*b,3*k)%1
            q=(-x-y)%1
            assert (k*q-Fraction(j,3)).denominator==1
            multisets.add(tuple(sorted((x,y,q))))
    delta=sum((k*a-j)%3==0 for a in range(3))
    expected=(k*k+3*k+2*delta)//6
    assert 6*expected==k*k+3*k+2*delta
    assert len(multisets)==expected
    all_equal=sum(len(set(q))==1 for q in multisets)
    two_equal=sum(len(set(q))==2 for q in multisets)
    distinct=sum(len(set(q))==3 for q in multisets)
    assert all_equal==delta
    assert two_equal==k-delta
    assert distinct==(k*k-3*k+2*delta)//6
    assert all_equal+3*two_equal+6*distinct==k*k
    return len(multisets)


def smooth_boundary_orbits(k):
    # Target roots q,q,q^-2 with q=exp(2*pi*i/7), a noncusp point.
    multisets=set()
    regular=0
    for a,b in combinations_with_replacement(range(k),2):
        x=Fraction(1+7*a,7*k)%1
        y=Fraction(1+7*b,7*k)%1
        q=(-x-y)%1
        assert (k*q+Fraction(2,7)).denominator==1
        multisets.add(tuple(sorted((x,y,q))))
        regular += a==b
    assert len(multisets)==k*(k+1)//2 and regular==k


def scaling_checks():
    # Check both the canonical real model and a nonreal rational scaling.
    p=lambda x:S.Poly(x,u,v,domain=S.QQ)
    rec=[p(3),p(u),p(u*u-2*v)]
    results=[]
    for k in range(2,21):
        if k>=3:
            rec.append(p(u)*rec[-1]-p(v)*rec[-2]+rec[-3])
        assert rec[k]==poly_T(k)
        assert poly_T(k).degree(v)==k//2
        results.append(k)
    d=S.Integer(1)+2*S.I
    dc=S.conjugate(d)
    beta=S.expand(-d*d/dc); alpha=S.conjugate(beta)
    pc=lambda x:S.Poly(x,u,v,domain=S.QQ_I)
    scaled=[pc(3),pc(d*u),pc(d*d*u*u+2*beta*dc*v)]
    for k in range(3,21):
        scaled.append(pc(d*u)*scaled[-1]+pc(beta*dc*v)*scaled[-2]
                      -pc(alpha*beta*beta)*scaled[-3])
    for k in range(2,21):
        assert scaled[k]==pc(d**k*poly_T(k).as_expr())
        target=-(-alpha/d)**k
        assert S.cancel(target**3+(alpha/beta)**(2*k))==0
    # Sharp radius: the largest possible discriminant at modulus r.
    r=S.symbols('r', real=True)
    assert S.expand((r-1)*(r+3)**3-(r**4+8*r**3+18*r*r-27))==0
    # Exact phase test for the earlier shifted model on the unit target circle.
    for k in range(2,31):
        for order in range(1,41):
            exceptional=(2*k)%order==0
            assert exceptional == (S.Rational(2*k,order).denominator==1)
    return results


def real_global_check(k,r):
    p=lambda x:S.Poly(x,z,w,t,domain=S.QQ)
    N=2*r+1
    Q=p(1)
    for h in range(1,r+1): Q*=p(z*w-h*h)
    f=p(z**N)+Q; fs=p(w**N)+Q
    B=p(-z**r); Bs=p(-w**r)
    A=p(w**(r+1))-B*p(z*w)
    G=A*f+B*fs
    C=p(t*S.Rational(1,k))*G
    D=p(t**3*S.Rational(1,k**3))*B*B*Bs
    L=[p(3),f,f*f+2*C]
    for h in range(3,k+1): L.append(f*L[-1]+C*L[-2]-D*L[-3])
    P=L[k]  # Crucially, no final conjugate-parameter power is added.
    n=k*N; m=k*r
    assert P.degree(z)==n and P.degree(w)<=m
    assert P.coeff_monomial(z**n)==1
    assert P.coeff_monomial(z**m*w**m)==1
    assert all(i+j<=n for (i,j,h),c in P.terms())
    if r>=2:
        assert all(i+j<n for (i,j,h),c in (P-p(z**n)).terms())
    else:
        assert all(h>0 for (i,j,h),c in P.terms() if i+j==n and (i,j)!=(n,0))
    return {'k':k,'r':r,'degree_z':n,'degree_w':int(P.degree(w)),
            'coefficient_field':'Q','pure_head':r>=2}


def real_circle_rank(m,n):
    # A real Laurent restriction A chosen coprime to its reciprocal.
    d=m+n
    x=S.symbols('x')
    A=S.Poly(x**d+2,x)
    revA=S.Poly(x**d*A.as_expr().subs(x,1/x),x)
    assert S.gcd(A,revA).degree()==0
    columns=[]
    for j in range(d+3):
        B=x**j
        # Numerator of 2i Im(Q bar P), after multiplying by x^(d+1).
        F=S.Poly(B*revA.as_expr()-x**(d+2-j)*A.as_expr(),x)
        columns.append([F.nth(d+1+h) for h in range(1,d+2)])
    matrix=S.Matrix.hstack(*(S.Matrix(col) for col in columns))
    assert matrix.rank()==d+1
    return [m,n,d+1]


def run():
    start=time.time()
    identities=[discriminant_checks(k) for k in range(2,16)]
    print('Discriminant-Jacobian identities for orders 2 through 15 passed.',flush=True)
    cusp_count=0
    for k in range(1,81):
        smooth_boundary_orbits(k)
        for j in range(3):
            cusp_orbits(k,j); cusp_count+=1
    scaling=scaling_checks()
    degrees=[real_global_check(k,r) for r in range(1,5) for k in range(2,9)]
    ranks=[real_circle_rank(m,n) for m in range(1,5) for n in range(m+1,m+5)]
    report={'status':'passed',
            'scope':'Finite exact diagnostics; all-degree fiber and existence theorems are proved in the TeX.',
            'discriminant_jacobian_identities':identities,
            'smooth_boundary_orders':80,'cusp_cases':cusp_count,
            'real_model_orders':scaling,
            'real_global_degree_cases':degrees,
            'real_circle_ranks':ranks,
            'seconds':round(time.time()-start,3)}
    path=Path(__file__).with_name('exact_local_models_verification.json')
    path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'passed','cusp_cases':cusp_count,
                      'global_cases':len(degrees),'circle_cases':len(ranks),
                      'seconds':report['seconds']}),flush=True)


if __name__=='__main__':
    run()
