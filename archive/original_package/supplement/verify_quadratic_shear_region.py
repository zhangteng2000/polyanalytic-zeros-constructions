"""Exact identities for the complete quadratic shear parameter region.
The fiber counts use the proved complete T_2 fiber theorem in the source.
"""
from pathlib import Path
import hashlib,json
import sympy as s
if not __debug__:raise RuntimeError('Run without -O: proof assertions must remain enabled')
x,y,u,v,r=s.symbols('x y u v r',real=True)
lam=x+s.I*y;bar=x-s.I*y
D=lambda a,b:s.expand(a*a*b*b+18*a*b-4*(a**3+b**3)-27)
c=lam**2+2*bar;cb=bar**2+2*lam
delta=s.expand((x*x+y*y)**2+18*(x*x+y*y)+4*(lam**3+bar**3)-27)
assert s.expand(D(c,cb)-(x*x+y*y-1)**2*delta)==0
assert s.expand(delta-D(-lam,-bar))==0
z=u+s.I*v;zb=u-s.I*v
q=-2*z-lam;qb=-2*zb-bar
assert s.expand(q*q-2*qb-c-4*(z*z+zb+lam*z))==0
assert s.expand(r**4+18*r*r+8*r**3-27-(r-1)*(r+3)**3)==0
for a in [s.S(1),s.S(-3)]:assert a*a+2*a==3
omega=(-1+s.I*s.sqrt(3))/2
for j in range(3):
    for a in [s.S(1),s.S(-3)]:
        l=s.simplify(a*omega**j)
        assert s.simplify(l*l+2*s.conjugate(l)-3*omega**(2*j))==0
report={'all_passed':True,'arithmetic':'exact symbolic identities',
        'four_simple_zeros':'Delta_b(eta)<0 and abs(eta)!=1',
        'two_simple_zeros':'Delta_b(eta)>0',
        'boundary':'3 distinct (2 regular), except (b/conj(b))*eta^3 in {1,-27}: 2 distinct (1 regular)',
        'unique_largest_disk_for_each_b':'center 0, radius 1',
        'uniform_four_zero_region_over_all_b':'abs(eta)<1',
        'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
dest=Path(__file__).with_name('quadratic_shear_region_verification.json')
dest.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2))
