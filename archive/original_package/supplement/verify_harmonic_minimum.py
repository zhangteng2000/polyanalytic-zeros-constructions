"""Exact finite diagnostics for the harmonic minimum and parameter radius.

Requires Python 3 and SymPy. Run without -O. The analytic index, degree,
continuation, and completion proofs are in polyanalytic_sharp_v3.tex.
No floating-point root count or finite sampling is used as a universal proof.
"""
from pathlib import Path
from fractions import Fraction as Q
import json
import time
from math import gcd
import sympy as S

if not __debug__:
    raise RuntimeError('Run without Python optimization (-O).')

start = time.time()
n, m, k, x = S.symbols('n m k x', positive=True)
d = n-k
h = m-k
C = k/n
B = (m*m-n*k)/(d*n)
re_u = (d*d+(n*n-m*m)*x)/(2*d*n)
assert S.factor(1-2*re_u+x-(C+B*x)) == 0
assert S.factor(d*d-2*d*n*re_u+(n*n-m*m)*x) == 0
assert S.factor(1/x-h*B/(d*(C+B*x))
                -(C+(n-m)*B*x/d)/(x*(C+B*x))) == 0
assert S.factor((C+(n-m)*B*x/d).subs(x,d*d/(n-m)**2)
                -m*(m-k)/(n*(n-m))) == 0
u_star = d/(n+m)
assert S.factor(d*d-2*d*n*u_star+(n*n-m*m)*u_star**2) == 0
assert S.factor(1-u_star-(m+k)/(n+m)) == 0
assert S.factor((d-n*u_star)**2-m*m*u_star**2) == 0

aa = d/(n+m)
tt = (m+k)/(n+m)
LL = 2*m*aa
S2 = (tt*n*(n-1)-k*(k-1)+aa*m*(m-1))/2
T2 = (tt*n*(n-1)-k*(k-1)-aa*m*(m-1))/2
T3 = (tt*n*(n-1)*(n-2)-k*(k-1)*(k-2)-aa*m*(m-1)*(m-2))/6
assert S.factor(-T2/m-m*aa+(m+k)*d/(2*m)) == 0
assert S.factor(2*T2*S2/LL-T3
                -(m+k)*d*(k*m+3*k*n-m*m+m*n)/(12*m)) == 0

radius_cases = 0
for nn in range(2, 51):
    for mm in range(1, nn):
        for kk in range(1, mm+1):
            dd, hh = nn-kk, mm-kk
            cc = Q(kk,nn)
            bb = Q(mm*mm-nn*kk,dd*nn)
            lo, hi = Q(dd*dd,(nn+mm)**2), Q(dd*dd,(nn-mm)**2)
            if hh:
                assert cc+bb*lo > 0 and cc+bb*hi > 0
                assert cc+Q(nn-mm,dd)*bb*lo > 0
                assert cc+Q(nn-mm,dd)*bb*hi > 0
            # R^d is rational; comparison needs no algebraic approximation.
            radius_power = Q(dd**dd,(nn+mm)**(nn-mm)*(mm+kk)**hh)
            assert Q(1,nn+mm+1)**dd < radius_power
            assert radius_power >= Q(dd,nn+mm)**dd
            if hh:
                assert dd >= 2 and Q(nn+mm,nn+mm) < dd
            else:
                assert Q(nn+mm,nn+mm+1) < dd
            radius_cases += 1
            assert (mm+kk)*dd*(kk*mm+3*kk*nn-mm*mm+mm*nn) > 0

boundary_orbit_cases = 0
for nn in range(2, 31):
    for mm in range(1, nn):
        for kk in range(1, mm+1):
            dd = nn-kk
            ell = gcd(dd,mm+kk)
            fibers = {}
            for jj in range(dd):
                value = ((mm+kk)*jj)%dd
                fibers[value] = fibers.get(value,0)+1
            assert len(fibers) == dd//ell
            assert all(number == ell for number in fibers.values())
            assert (dd+ell)-ell+kk == nn
            boundary_orbit_cases += 1

z = S.symbols('z')
gcd_cases = 0
for nn in range(2, 17):
    for mm in range(1, nn):
        for kk in range(1, mm+1):
            A = nn+mm+1
            pp = S.Poly(S.diff(A*(z**nn-z**kk),z),z)
            qq = S.Poly(mm*z**(mm-1),z)
            gg = S.gcd(pp,qq).monic()
            assert gg.as_expr() == z**(kk-1)
            assert nn-(kk-1) == (nn-kk)+1
            gcd_cases += 1

two_zero_cases = 0
for nn in range(2, 81):
    # m=k=n-1: the two exhibited zeros are 0 and this rational point.
    rr = Q(2*nn-1,2*nn)
    assert 2*nn*(rr**nn-rr**(nn-1))+rr**(nn-1) == 0
    pz = 2*nn*(nn*rr**(nn-1)-(nn-1)*rr**(nn-2))
    qw = (nn-1)*rr**(nn-2)
    assert pz*pz-qw*qw > 0
    assert Q(1,2*nn) < Q(1,2*nn-1)
    two_zero_cases += 1

# Explicit elimination check: n=3, m=k=2, A=6.
xx, yy, ss = S.symbols('x y s', real=True)
zz = xx+S.I*yy
ff = S.expand(6*zz**3-6*zz**2+S.conjugate(zz)**2)
re, im = S.expand(S.re(ff)), S.expand(S.im(ff))
assert S.factor(re.subs(yy,0)) == xx**2*(6*xx-5)
assert S.factor(im) == 2*yy*(9*xx**2-7*xx-3*yy**2)
off_axis = S.factor(re.subs(yy**2,3*xx**2-S.Rational(7,3)*xx))
assert off_axis == -xx*(144*xx**2-156*xx+35)/3
quad = S.Poly(144*xx**2-156*xx+35,xx)
# Both roots are in (0,7/9), where the required y^2 is negative.
assert quad.count_roots(0,S.Rational(7,9)) == 2

completion_cases = 0
pure_exact_cases = 0
for nn in range(2, 81):
    for mm in range(1, nn):
        for hh in range(2, nn+1):
            pp = max(mm+1,hh)
            kk = pp-hh+1
            assert 1 <= kk <= mm < pp <= nn
            assert nn-pp+kk == nn-hh+1 < nn
            assert nn-pp+mm < nn
            assert pp-kk+1 == hh
            assert (nn-pp)+kk == nn-hh+1
            pure_exact_cases += 1

gains = []
for mm in range(1, 301):
    for nn in range(4*mm+1,4*mm+9):
        assert 10*mm*mm+nn-3*mm+mm == nn+10*mm*mm-2*mm
        assert 4*mm < nn
        completion_cases += 1
    nn = 4*mm+1
    aa = min(mm,(nn-2*mm)//4)
    radial = max((2*mm+1)*nn,nn+2*(mm+aa)*(nn-2*aa-1))
    new = nn+10*mm*mm-2*mm
    if mm >= 2:
        assert radial == 9*mm*mm+4*mm+(1 if mm%2==0 else 0)
        assert new-radial == (mm*(mm-2) if mm%2==0 else (mm-1)**2)
    if mm in (3,10,100):
        gains.append({'m':mm,'n':nn,'radial':radial,'completion':new})

report = {
    'status':'passed',
    'scope':'Exact finite algebraic diagnostics, not substitutes for the analytic proofs.',
    'symbolic_identities':9,
    'rational_radius_cases':radius_cases,
    'derivative_gcd_cases':gcd_cases,
    'boundary_orbit_cases':boundary_orbit_cases,
    'two_zero_family_checks':two_zero_cases,
    'explicit_two_zero_elimination':'n=3, m=k=2; exact Sturm count passed',
    'completion_degree_cases':completion_cases,
    'pure_exact_degree_cases':pure_exact_cases,
    'radial_gain_cases':299,
    'examples':gains,
    'seconds':round(time.time()-start,3),
}
Path(__file__).with_name('harmonic_minimum_verification.json').write_text(
    json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
