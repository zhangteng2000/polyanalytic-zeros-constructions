"""Exact checks for the paper; these finite checks are not universal proofs.

Run from any directory with Python >=3.10 and SymPy >=1.13.
Reports are written beside this script. No archived certificate is required.
"""
from pathlib import Path
import argparse
import json
import math
import platform
import time
import sympy as sp

OUT = Path(__file__).resolve().parent
u = sp.symbols("u", positive=True)
r = sp.symbols("r", real=True)
records = []
start = time.monotonic()

def record(name, **data):
    records.append({"check": name, **data})
    print(name, json.dumps(data, ensure_ascii=True), flush=True)

def positive_count(expr, variable=r, start=0):
    p = sp.Poly(sp.expand(expr), variable, extension=True)
    count = int(p.count_roots(start, sp.oo))
    assert sp.gcd(p, p.diff()).degree() == 0, "Repeated algebraic root"
    return count

def check_m1():
    # u=e^t>1. Clearing the positive denominator gives the exact
    # selected-branch equation at epsilon=1/n^2.
    for n in (2, 3, 4, 5, 6, 8, 10):
        counts = []
        for j in range(1, n):
            sine = sp.sin(sp.pi*j/n)
            polynomial = u**(2*n) + 1 - 2*n*n*sine*(u**(n+1)-u**(n-1))
            counts.append(positive_count(polynomial, u, 1))
        assert counts == [2]*(n-1)
        record("m1_exact_Sturm", n=n, branch_counts=counts,
               total=n+sum(counts), expected=3*n-2)

def check_general_chebyshev():
    # On z_j(t), express z in u=e^t. The imaginary part of T_m(z-a)
    # is a Laurent polynomial with real algebraic coefficients.
    epsilon = sp.Rational(1, 10**10)
    a = 3
    for n, m in ((4, 2), (4, 3), (5, 2), (6, 2), (6, 3), (6, 5)):
        upper, lower = [], []
        for j in range(1,n):
            c, s = sp.cos(sp.pi*j/n), sp.sin(sp.pi*j/n)
            z = c*(u+1/u)/2 + sp.I*s*(u-1/u)/2
            translated = sp.expand(sp.chebyshevt(m,z-a))
            imaginary = sp.expand((translated-sp.conjugate(translated))/(2*sp.I))
            for sign, counts in ((1,upper),(-1,lower)):
                equation = epsilon*(-1)**j*(u**n+u**(-n))/2-sign*2*imaginary
                cleared = sp.expand(equation*u**n)
                counts.append(positive_count(cleared,u,1))
        expected = (m+2)*n-m-math.gcd(m,n)
        assert n+sum(upper)+sum(lower)==expected
        record("Chebyshev_exact_Sturm", n=n,m=m,a=a,
               epsilon=str(epsilon),upper=upper,lower=lower,total=expected)

def check_radial():
    for a,b,k in ((0,0,1),(0,1,3),(0,2,5),(1,1,3),
                  (2,1,4),(2,2,5),(3,2,6),(3,3,7)):
        eta=sp.Rational(1,8*(a+b+2))
        B=sum((-1)**j*eta**(4*j*j)*r**(2*j) for j in range(b+1))
        A=sum((-1)**(b+1+j)*eta**((k+2*j)**2)*r**(2*j) for j in range(a+1))
        plus=positive_count(B+r**k*A)
        minus=positive_count(B-r**k*A)
        assert plus==a+b+1 and minus==a+b
        assert sp.gcd(sp.Poly(A,r),sp.Poly(B,r)).degree()==0
        record("radial_exact_Sturm",a=a,b=b,k=k,eta=str(eta),
               plus=plus,minus=minus,planar_count=k*(plus+minus))

def check_quartic():
    plus=r**3+r**2-sp.Rational(1,8)
    minus=-r**3+r**2-sp.Rational(1,8)
    radii=[(sp.sqrt(5)-1)/4,sp.Rational(1,2),(1+sp.sqrt(5))/4]
    for radius,poly in zip(radii,[plus,minus,minus]):
        assert sp.simplify(poly.subs(r,radius))==0
        assert sp.simplify(sp.diff(poly,r).subs(r,radius))!=0
    assert positive_count(plus)==1 and positive_count(minus)==2
    record("quartic_exact",radii=[str(x) for x in radii],origin_jacobian="-1/64",total=16)

def sign_exact(x):
    y=sp.simplify(x)
    if y.is_positive: return 1
    if y.is_negative: return -1
    raise AssertionError("Unable to establish exact sign: "+str(y))

def variations(signs):
    return sum(a!=b for a,b in zip(signs,signs[1:]))

def check_harmonic_angles():
    # Exact signs of real algebraic numbers, no floating-point decisions.
    for n in range(2,11):
        for m in range(1,n):
            g=math.gcd(n,m)
            delta=g if (n//g)%2 and (m//g)%2 else 0
            sigma=-1 if delta==g and n//g%4==1 else 1
            exponents=list(range(m+1))+[n]
            count=0
            for j in range(n):
                signs=[sigma*sign_exact(sp.cos(sp.pi*k*j/n)
                       -sp.Rational(1,4*n)*sp.sin(sp.pi*k*j/n))
                       for k in range(m+1)]+[(-1)**j]
                count+=variations(signs)
                count+=variations([s*(-1)**e for s,e in zip(signs,exponents)])
            assert count==n*(m+1)+delta,(m,n,count)
            record("restricted_harmonic_exact_signs",n=n,m=m,total=count,delta=delta)

def check_formulas():
    for m in range(1,101):
        for n in range(m+1,501):
            bound=min(n*n+m*m,(4*m+1)*n-5*m)
            if n<=2*m+1:
                lower=n*n+m*m
            else:
                s=max(0,3*m+1-n)
                k=n-m
                lower=(k+s)*(4*m-2*s+1)+s*s
                assert lower==(4*m+1)*(n-m)+max(0,3*m-n)*max(0,3*m+1-n)
            assert 0<=bound-lower<=4*m*(m-1)
            assert (bound-lower)%2==0
    record("integer_formula_grid",m_max=100,n_max=500,status="pass")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--group",choices=("all","m1","chebyshev","radial","angles"),default="all")
    args=parser.parse_args()
    if args.group in ("all","m1"):check_m1()
    if args.group in ("all","chebyshev"):check_general_chebyshev()
    if args.group in ("all","radial"):
        check_radial()
        check_quartic()
        check_formulas()
    if args.group in ("all","angles"):check_harmonic_angles()
    report={"status":"PASS","python":platform.python_version(),"sympy":sp.__version__,
            "arithmetic":"exact rational and real algebraic arithmetic",
            "scope":"Finite consistency checks; universal claims rely on the manuscript proofs.",
            "seconds":round(time.monotonic()-start,3),"records":records}
    name="verification_report.json" if args.group=="all" else "verification_"+args.group+".json"
    (OUT/name).write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("PASS:",len(records),"checks",flush=True)

if __name__=="__main__":main()
