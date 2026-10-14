"""Exact algebra supporting the analytic sharp constants for powers 2--5.

No numerical maximizer or root finder is used for any pass/fail decision.
The inequalities reducing the compact trace domain are proved in the note.
Run with Python and SymPy; write the report beside this file.
"""
import sys
if not __debug__:
    raise RuntimeError('Do not run this exact verifier with Python -O.')
from pathlib import Path
import hashlib,json
import sympy as s

u,v,r,t,q=s.symbols('u v r t q')
checks=[]
def zero(name,expr):
    assert s.cancel(s.expand(expr))==0,name
    checks.append(name)
def yes(name,condition):
    assert condition,name
    checks.append(name)
T=[s.Integer(3),u,u*u-2*v]
for k in range(3,6):T.append(s.expand(u*T[-1]-v*T[-2]+T[-3]))
zero('second power identity',T[2]-(u*u-2*v))
zero('third power divided identity',T[3]-3-u*(u*u-3*v))
zero('fourth power identity',T[4]-(u**4-4*u*u*v+2*v*v+4*u))
zero('fifth power identity',T[5]-(u**5-5*u**3*v+5*u*v*v+5*u*u-5*v))
e1,e2,e3=s.symbols('e1 e2 e3')
zero('third power pair difference identity',(e1+e2+e3)**2-3*(e1*e2+e1*e3+e2*e3)-((e1-e2)**2+(e1-e3)**2+(e2-e3)**2)/2)
zero('unit pair distance identity',sum((a-b)*(1/a-1/b) for a,b in [(e1,e2),(e1,e3),(e2,e3)])-(9-(e1+e2+e3)*(1/e1+1/e2+1/e3)))
zero('second power polar expression',T[2].subs({u:r*q,v:r/q})/(r/q)-(r*q**3-2))
zero('third power polar expression',(T[3].subs({u:r*q,v:r/q})-3)/(r*q)-(r*r*q*q-3*r/q))
four=4-4*r*r+r**3*q**3+2*r/q**3
zero('fourth power polar expression',T[4].subs({u:r*q,v:r/q})/(r*q)-four)
f4=4+2*r-4*r*r+r**3
rho4=(4-s.sqrt(10))/3
C4=(52+20*s.sqrt(10))/27
zero('fourth power critical point',s.diff(f4,r).subs(r,rho4))
zero('fourth power exact maximum',f4.subs(r,rho4)-C4)
yes('fourth power unique derivative sign change',0<rho4<1 and (4+s.sqrt(10))/3>1)
yes('fourth power strict improvement over 4',C4>4)
a=-5+5*r*r;b=5*r*(1-r*r);c=r**4
zero('fifth power polar expression',q**2*T[5].subs({u:r*q,v:r/q})/(r*q)-(a+b*q**3+c*q**6))
square=(a-c)**2+b*b+2*b*(a+c)*t+4*a*c*t*t
zero('quadratic modulus identity',(a+b*q+c*q*q)*(a+b/q+c/q**2)-square.subs(t,(q+1/q)/2))
vertex=-(5-5*r*r-r**4)/(4*r**3)
zero('fifth power angular vertex',s.diff(square,t).subs(t,vertex))
K=5-5*r*r-4*r**3-r**4
zero('vertex threshold',vertex+1+K/(4*r**3))
yes('vertex endpoint positivity',K.subs(r,s.Rational(3,5))==s.Rational(1379,625))
zero('vertex monotonicity',s.diff(K,r)-(-10*r-12*r*r-4*r**3))
f5=5+5*r-5*r*r-5*r**3-r**4
zero('fifth power maximum at negative phase',f5-(-a+b-c))
g=4*r**3+15*r*r+10*r-5
zero('fifth power radial derivative',s.diff(f5,r)+g)
zero('fifth power strict concavity',s.diff(f5,r,2)-(-10-30*r-12*r*r))
yes('fifth power root isolation signs',g.subs(r,s.Rational(3,10))<0<g.subs(r,s.Rational(1,3)))
yes('fifth power exceeds outer-domain bound',f5.subs(r,s.Rational(3,10))>5)
zero('fifth power equality trace polynomial',g-s.Poly(g,r).as_expr())
report={'all_exact_checks_passed':True,'exact_check_count':len(checks),'checks':checks,
 'C2':'3','C3':'4','C4':'(52+20*sqrt(10))/27',
 'C5':'5+5*rho-5*rho^2-5*rho^3-rho^4; 4*rho^3+15*rho^2+10*rho-5=0; 3/10<rho<1/3',
 'scope':'Exact algebra supports the analytic proof, not an independent theorem prover; no numerical search certifies an endpoint.',
 'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
destination=Path(__file__).with_name('low_order_power_constants_verification.json')
destination.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2))
