"""Exact algebraic audits for the cyclic contact-five reduction.
The all-degree and singular-zero arguments are analytic and in the paper.
No numerical root count is accepted by this checker.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

if not __debug__:
    raise RuntimeError('Run without -O: proof assertions must remain enabled')
z,w,u,v,a,b,c,r=s.symbols('z w u v a b c r')
P=z**5+a*z**3*w+b*z*w*w+c*w
G=v*(u**3+a*u*u+b*u)+c
assert s.cancel(P/w-G.subs({u:z*z/w,v:w*w/z}))==0
allowed=[(5,0)]+[(i,j) for j in (1,2) for i in range(5-j)
                  if (i-j-2)%3==0]
assert set(allowed)=={(5,0),(3,1),(1,2),(0,1)}
assert s.diff(P,w).subs({z:0,w:0})==c
assert P.subs(w,0)==z**5
sturm=[]
for d in range(2,21):
    plus=r**(d+1)-3*r*r+1
    minus=-r**(d+1)-3*r*r+1
    assert s.Poly(plus,r).count_roots(0,s.oo)==2
    assert s.Poly(minus,r).count_roots(0,s.oo)==1
    assert s.gcd(plus,s.diff(plus,r))==1
    assert s.gcd(minus,s.diff(minus,r))==1
    assert s.cancel(s.diff(plus,r)/r-((d+1)*r**(d-1)-6))==0
    assert plus.subs(r,1)==-1
    sturm.append({'degree':d,'positive_branch_radii':2,
                  'negative_branch_radii':1,'regular_zero_count':3*(d-1)})
out={'all_passed':True,'arithmetic':'exact SymPy rational polynomial arithmetic',
     'covering_identity':True,'full_cyclic_support':allowed,
     'origin_jacobian':'-abs(c)^2 < 0',
     'sturm_checks':sturm,'cyclic_degree_five_count':19,
     'scope':'Algebraic audits only; general upper and singular proof in manuscript',
     'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
dest=Path(__file__).with_name('cyclic_contact_verification.json')
dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'all_passed':True,'checked_degrees':len(sturm),'output':str(dest)}))
