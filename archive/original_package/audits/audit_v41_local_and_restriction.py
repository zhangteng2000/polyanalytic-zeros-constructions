"""Exact diagnostics for the analytic selective-splitting proof, not a root search."""
from pathlib import Path
import sympy as s,json
x,y,l,z=s.symbols('x y lambda z',real=True)
u=x+s.I*y;Q=s.expand(u**3+l*u**2*s.conjugate(u)-l*u+s.conjugate(u))
R=s.re(Q).expand();I=s.im(Q).expand()
assert s.expand(R-x*((1+l)*x*x+(l-3)*y*y-l+1))==0
assert s.expand(I-y*((3+l)*x*x+(l-1)*y*y-l-1))==0
J=s.expand(s.diff(R,x)*s.diff(I,y)-s.diff(R,y)*s.diff(I,x))
def at(X,Y):
    return s.factor(s.Poly(J,x,y).eval({x:s.sqrt(X),y:s.sqrt(Y)}))
assert at(0,0)==l*l-1 or s.expand(at(0,0)-(l*l-1))==0
assert s.factor(at((l-1)/(l+1),0)+8*(l-1)/(l+1))==0
assert s.factor(at(0,(l+1)/(l-1))+8*(l+1)/(l-1))==0
assert at(s.Rational(1,2),s.Rational(1,2))==8
assert s.Matrix([[l+1,l-3],[l+3,l-1]]).det()==8
rank_results=[]
for d in range(3,10):
    for kind in ['rectangle_even','rectangle_odd','pure_even','pure_odd']:
        if kind.startswith('rectangle'):
            p=z**(d-1)+z+1;q=z**d+2
            n=2*d+(kind.endswith('odd'));cap=n+2*d
            terms=[z**i*(-p)**j*q**(2-j) for j in range(3) for i in range(n+1)]
        else:
            p=z**d+z+1;q=z**(d-2)+2
            D=2*d-(kind.endswith('even'));cap=D+2*d-2
            terms=[z**i*(-p)**j*q**(2-j) for j in range(3) for i in range(D-j+1)]
        assert s.gcd(p,q)==1
        polys=[s.Poly(t,z) for t in terms]
        M=s.Matrix([[p.nth(i) for p in polys] for i in range(cap+1)])
        rank=M.rank();assert rank==cap+1
        rank_results.append(dict(d=d,kind=kind,rank=rank))
for d in range(4,1001):
    N=5*d-5
    se=min(N,4*d);so=min(N-1,4*d)
    assert 2*N+2*se<=18*d-10
    if d>=5:assert 2*N+2*se==18*d-10
    if d>=6:assert 2*N+2*so+7==18*d-3
    N=5*d-6
    assert 4*d-3<=N-1
    assert 2*N+2*(4*d-3)==9*(2*d)-18
    assert 2*N+2*(4*d-3)+7==9*(2*d+1)-20
out={'all_passed':True,'nine_zero_model':'exact real/imaginary equations and all four Jacobian formulas',
     'restriction_rank_examples':rank_results,'degree_and_count_checks':'d=4,...,1000',
     'scope':'Finite exact diagnostics only; all-degree surjectivity, zero confinement, and global counts are proved analytically in the manuscript.'}
path=Path(__file__).with_name('v41_local_and_restriction_audit.json')
path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='restriction_rank_examples'},indent=2))
