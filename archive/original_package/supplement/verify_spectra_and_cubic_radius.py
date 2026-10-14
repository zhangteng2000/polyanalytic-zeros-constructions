"""Exact diagnostics for the spectrum and optimal cubic-radius proofs.

Requires Python and SymPy, run without -O. The general transversality and
continuation arguments are analytic proofs in the TeX, not finite tests.
"""
from pathlib import Path
from fractions import Fraction as F
import json
import time
import sympy as s

if not __debug__:
    raise RuntimeError('Run without -O: this verifier relies on assertions.')


def cubic_identities():
    x,k,r=s.symbols('x k r',positive=True)
    c=(k*k-x*x-8*x**3)/(2*s.sqrt(x)**3*(3*x+k))
    psi=k*k+(x-k)**2*(2-k-5*x)/(k+3*x)
    e2=k*k/x-2*k+x+x*x-2*s.sqrt(x)*(k-x)*c
    assert s.simplify(e2-psi)==0
    A=4*x**3-x*x-2*k*x-k*k
    B=16*x**3-x*x+2*k*x-k*k
    assert s.expand((k*k-x*x-8*x**3)**2-4*x**3*(3*x+k)**2-A*B)==0
    assert s.expand(B-4*A-(3*x*x+10*k*x+3*k*k))==0
    derivative=-2*(x-k)*(15*x*x+(9*k-3)*x-5*k)/(k+3*x)**2
    assert s.factor(s.diff(psi,x)-derivative)==0
    kr=r*r*(2*r-1)
    assert s.factor(A.subs({x:r*r,k:kr}))==0
    assert s.factor(psi.subs({x:r*r,k:kr})-(2*r-3*r*r)**2)==0
    assert s.factor(kr-(2*r-3*r*r))==2*r*(r*r+r-1)
    optimum=(5*s.sqrt(5)-11)/2
    assert s.simplify(optimum**2+11*optimum-1)==0
    assert s.simplify(kr.subs(r,(s.sqrt(5)-1)/2)-optimum)==0
    assert optimum>s.Rational(9,100)
    assert optimum<s.Rational(91,1000)
    # Real singular witness v=-r and e=2r-3r^2.
    ev=2*r-3*r*r
    assert s.expand(-r**3+r*r-r*ev-kr)==0
    assert s.expand((3*r*r-r)**2-(-r+ev)**2)==0
    # Exact rational sampling of the all-x lower bound (diagnostic only).
    count=0
    for i in range(1,120):
        rr=F(1,2)+F(i,720)
        kk=rr*rr*(2*rr-1)
        bound=min(kk,2*rr-3*rr*rr)**2
        for j in range(121):
            xx=rr*rr*F(j,120)
            val=kk*kk+(xx-kk)**2*(2-kk-5*xx)/(kk+3*xx)
            assert val>=bound
            count+=1
    return {'rational_radius_inequalities':count,
            'optimal_kappa':str(optimum),'identities':'passed'}


def spectrum_algebra():
    a,b=s.symbols('a b')
    H=s.Matrix([[1,a,a*a,a**3],[0,1,2*a,3*a*a],
                [1,b,b*b,b**3],[0,1,2*b,3*b*b]])
    assert s.factor(H.det())==(a-b)**4
    x,y,z,w=s.symbols('x y z w',real=True)
    P2=z*z*w-s.Rational(5,2)*z-s.Rational(1,2)*w-2
    P4=z*z*w-4*z-w-2*s.I
    examples=[]
    for P,roots,fold in [(P2,[(2,0),(-1,0)],(-1,0)),
                         (P4,[(0,2),(0,-1),(2,1),(-2,1)],(0,-1))]:
        C=s.expand(P.subs({z:x+s.I*y,w:x-s.I*y}))
        U=s.re(C).expand(); V=s.im(C).expand()
        jac=s.Matrix([U,V]).jacobian([x,y])
        for px,py in roots:
            assert U.subs({x:px,y:py})==V.subs({x:px,y:py})==0
            mat=jac.subs({x:px,y:py})
            assert mat.rank()==(1 if (px,py)==fold else 2)
        if P==P2:
            assert s.factor(U.subs(y,0))==(x-2)*(x+1)**2
            assert s.factor(V)==y*(x*x+y*y-2)
            assert s.diff(U,x,2).subs({x:-1,y:0})==-6
        else:
            assert s.factor(V.subs(x,0))==(y-2)*(y+1)**2
            assert s.factor(U)==x*(x*x+y*y-5)
            assert s.diff(V,y,2).subs({x:0,y:-1})==-6
        examples.append({'count':len(roots),'roots':roots,'fold':fold})
    # Rational correction enforces a zero and Jacobian zero identically.
    z,w=s.symbols('z w')
    Q=z**4+z*z*w+(1+2*s.I)*z*w*w+3*w+2-s.I
    count=0
    for i in range(1,31):
        aa=s.Rational(i,7)+s.I*s.Rational(i+1,11)
        ac=s.conjugate(aa)
        t=s.Rational(i,13)
        gamma=(1-t*t+2*s.I*t)/(1+t*t)
        Qw=s.expand(s.diff(Q,w).subs({z:aa,w:ac}))
        Qz=s.expand(s.diff(Q,z).subs({z:aa,w:ac}))
        A=s.expand(gamma*Qw-Qz)
        B=s.expand(-Q.subs({z:aa,w:ac})-A*aa)
        R=s.expand(Q+A*z+B)
        assert s.expand(R.subs({z:aa,w:ac}))==0
        rz=s.expand(s.diff(R,z).subs({z:aa,w:ac}))
        rw=s.expand(s.diff(R,w).subs({z:aa,w:ac}))
        assert s.expand(rz*s.conjugate(rz)-rw*s.conjugate(rw))==0
        assert rw!=0
        count+=1
    return {'hermite_determinant':'(a-b)^4',
            'low_degree_fold_examples':examples,
            'rational_zero_and_jacobian_checks':count}


def run():
    start=time.time()
    report={'scope':'Exact finite diagnostics; analytic proofs remain in the TeX.',
            'cubic_radius':cubic_identities(),'spectra':spectrum_algebra(),
            'status':'passed','seconds':round(time.time()-start,3)}
    Path(__file__).with_name('spectra_and_cubic_radius_verification.json').write_text(
        json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    run()
