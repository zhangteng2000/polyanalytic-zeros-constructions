"""Exact symbolic identities for the optimal cubic shear disk.

This checks the algebra used by the all-parameter analytic proof in
polyanalytic_sharp_v3.tex. It is not a numerical root-counting argument.
Requires SymPy. Run without Python's -O option.
"""
import argparse
import json
from pathlib import Path
import sympy as s


def main():
    if not __debug__:
        raise RuntimeError('Run without -O: exact checks must remain enabled.')
    x, c, t, a = s.symbols('x c t a', positive=True)
    z, w, e = s.symbols('z w e')
    checks = []

    def check(name, lhs, rhs=0):
        assert s.simplify(lhs-rhs) == 0, name
        checks.append(name)

    F = z**3-3*z*w-3*e*z**2+c
    E = (z**3-3*z*w+c)/(3*z**2)
    check('parameter elimination', F.subs(e, E))
    check('critical differential', s.diff(F,z).subs(e,E), z**2+3*w-2*c/z)
    cosine = -(x**3-12*c*x+4*c**2)/(2*x**s.Rational(3,2)*(3*x-2*c))
    crit = x**3+(3*x-2*c)**2+2*x**s.Rational(3,2)*(3*x-2*c)*cosine-9*x**2
    check('critical-angle elimination',crit)
    Q = (-2*c**3+13*c**2*x-c*x**3-24*c*x**2+2*x**4+9*x**3)/(3*x**2*(3*x-2*c))
    modulus = x/9+(c/(3*x)-1)**2+2*s.sqrt(x)/3*(c/(3*x)-1)*cosine
    check('critical parameter modulus', Q, modulus)
    A = x**4-c*x**3+9*c*x**2-10*c**2*x+4*c**3
    check('critical-modulus derivative',s.diff(Q,x),2*(3*x-c)*A/(3*x**3*(3*x-2*c)**2))
    B = x**3-36*x**2+24*c*x-4*c**2
    check('angular feasibility factorization',
          4*x**3*(3*x-2*c)**2-(x**3-12*c*x+4*c**2)**2,
          -(x**3-4*c**2)*B)
    check('positive quartic scaling', A.subs(x,c*t)/c**3,
          c*t**3*(t-1)+9*t**2-10*t+4)
    check('positive quartic certificate for 0<t<1', A.subs(x,c*t)/c**3,
          (2*t**2-t)**2+8*(t-s.Rational(5,8))**2+s.Rational(7,8)
          +(4-c)*t**3*(1-t))
    check('positive quadratic for t>=1',9*t**2-10*t+4,
          9*(t-s.Rational(5,9))**2+s.Rational(11,9))
    check('feasibility at c/3',B.subs(x,c/3),c**3/27)
    check('feasibility at 2c/3',B.subs(x,2*c/3),4*c**2*(2*c/27-1))
    check('outer endpoint factor',B.subs({c:4*a**3,x:4*a**2}),
          12*(4*a**2)*(8*a**3-12*a**2))
    check('inner minimum',Q.subs(x,c/3),c/27)
    check('outer minimum',Q.subs({c:4*a**3,x:4*a**2}),(1-a)**2)
    check('inner critical angle',cosine.subs(x,c/3),s.sqrt(c)/(6*s.sqrt(3)))
    check('outer critical angle',cosine.subs({c:4*a**3,x:4*a**2}),1)
    check('base cubic discriminant',(3-c)**4+18*(3-c)**2-8*(3-c)**3-27,-c**3*(4-c))
    check('optimal inner radius squared',(c/27).subs(c,s.Rational(27,16)),s.Rational(1,16))
    check('optimal outer radius',1-a.subs(a,s.Rational(3,4)),s.Rational(1,4))

    # Real Taylor coordinates at an inner nearest singularity.
    X,Y,r,b=s.symbols('X Y r b',real=True)
    U=r+X+s.I*Y
    phase=r/6-s.I*b
    shifted=s.expand(U**3-r*U**2+phase*(3*r**2-3*((r+X)**2+Y**2)))
    real=s.re(shifted).expand()
    imag=s.im(shifted).expand()
    check('inner fold first x component',s.diff(real,X).subs({X:0,Y:0}))
    check('inner fold first y component',s.diff(real,Y).subs({X:0,Y:0}))
    check('inner rank-one image x derivative',s.diff(imag,X).subs({X:0,Y:0}),6*r*b)
    check('inner rank-one image y derivative',s.diff(imag,Y).subs({X:0,Y:0}),r**2)
    q2=sum(co*X**i*Y**j for (i,j),co in s.Poly(real,X,Y).terms() if i+j==2)
    check('inner fold real quadratic',q2,r*(3*X**2-5*Y**2)/2)
    check('inner fold kernel coefficient',s.expand(q2.subs({X:r,Y:-6*b})).subs(b**2,1-r**2/36),
          4*r**3-90*r)
    # Outer singularity: v=t and c=t^3/2, eta=t/2-1.
    V=t+X+s.I*Y
    outer=s.expand(V**3-3*((t+X)**2+Y**2)-3*(t/2-1)*V**2+t**3/2)
    check('outer root',outer.subs({X:0,Y:0}))
    check('outer x derivative',s.diff(s.re(outer),X).subs({X:0,Y:0}))
    check('outer y derivative',s.diff(s.im(outer),Y).subs({X:0,Y:0}),6*t)
    check('outer kernel second derivative',s.diff(s.re(outer),X,2).subs({X:0,Y:0}),3*t)
    report={
        'status':'all exact symbolic identities passed',
        'checks':checks,
        'radius':'min(sqrt(c/27), 1-(c/4)^(1/3)), 0<c<4',
        'unique_optimal_c':'27/16',
        'optimal_radius':'1/4',
        'canonical_c':'3',
        'canonical_radius':'1-(3/4)^(1/3)',
        'boundary':'six inner folds, three outer folds, or all nine parameter values at the optimum',
        'scope':'Algebraic verification of an analytic all-parameter proof; no numerical search.'
    }
    parser=argparse.ArgumentParser()
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    if args.report:
        args.report.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
