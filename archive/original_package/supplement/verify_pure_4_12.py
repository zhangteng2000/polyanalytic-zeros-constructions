"""Verify 144 simple zeros using rational arithmetic only.

Run: python verify_pure_4_12.py [certificate.json] [--report report.json]

At each rational center c, let Y be the exact inverse real Jacobian.
On a square of radius r, the complex second derivatives give
||JF(z)-JF(c)||_infinity <= 4*r*(Mzz+Mzw+Mww).
Thus q=4*r*||Y||_infinity*(Mzz+Mzw+Mww) bounds the derivative
of z-Y*F(z). The exact inequality beta+q*r<r, where
beta=||Y*F(c)||_infinity, is a contraction certificate.
The pure highest part z^12 implies, with its conjugate w^12, no
intersection at infinity and no common projective component. Bezout's
144 intersections are exhausted by the 144 disjoint certified boxes.
Floating point is used only for optional human-readable summaries.
"""
from pathlib import Path
from fractions import Fraction as F
from math import lcm
import argparse
import json
import time
from functools import lru_cache

if not __debug__:
    raise RuntimeError('Run without Python optimization; exact assertions must be enabled.')

def p_add(*ps):
    out={}
    for p in ps:
        for mon,c in p.items():out[mon]=out.get(mon,F(0))+c
    return {mon:c for mon,c in out.items() if c}

def p_scale(c,p):return {mon:c*v for mon,v in p.items() if c*v}

def p_mul(a,b):
    out={}
    for (i,j),c in a.items():
        for (ii,jj),d in b.items():
            mon=(i+ii,j+jj); out[mon]=out.get(mon,F(0))+c*d
    return {mon:c for mon,c in out.items() if c}

def p_power(p,n):
    out={(0,0):F(1)}
    for _ in range(n):out=p_mul(out,p)
    return out

def reconstruct(kappa,tau):
    """Independent standard-library implementation of the order-four formula."""
    @lru_cache(None)
    def reduced_monomial(i,j):
        if i>=3:
            return p_add(p_scale(kappa,reduced_monomial(i-3,j)),
                         p_scale(F(-1),reduced_monomial(i-2,j+1)))
        if j>=3:
            return p_add(p_scale(kappa,reduced_monomial(i,j-3)),
                         p_scale(F(-1),reduced_monomial(i+1,j-2)))
        return {(i,j):F(1)}
    def reduction(p):
        return p_add(*(p_scale(c,reduced_monomial(i,j)) for (i,j),c in p.items()))
    f={(3,0):F(1),(1,1):F(1),(0,0):-kappa}
    g={(0,3):F(1),(1,1):F(1),(0,0):-kappa}
    A={(0,2):F(1)}; B={(1,0):F(-1)}
    D={(2,1):F(-1)}; H={(3,0):F(-1)}; C=D
    G=p_add(p_mul(A,f),p_mul(B,g))
    # Q4 = U4+U2*T+(H/4)U2+(A/4)U*T+((A*H-C)/16)U
    #      +T2/8+((A2+H)/16)T+A4/256+A2*H/64-A*C/64+H2/128.
    constant=reduction(p_add(p_scale(F(1,256),p_power(A,4)),
        p_scale(F(1,64),p_mul(p_power(A,2),H)),
        p_scale(F(-1,64),p_mul(A,C)),p_scale(F(1,128),p_power(H,2))))
    return p_add(p_power(f,4),p_scale(tau,p_mul(p_power(f,2),G)),
        p_scale(tau**2/4,p_mul(reduction(H),p_power(f,2))),
        p_scale(tau**2/4,p_mul(D,p_mul(f,G))),p_scale(tau**2/8,p_power(G,2)),
        p_scale(tau**3/16,p_mul(reduction(p_add(p_mul(A,H),p_scale(F(-1),C))),f)),
        p_scale(tau**3/16,p_mul(p_add(p_power(D,2),H),G)),
        p_scale(tau**4,constant))

def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def scale(c,a):return(c*a[0],c*a[1])

def evaluate(coeffs,x,y):
    z=(x,y); w=(x,-y)
    zp=[(F(1),F(0))]; wp=[(F(1),F(0))]
    for _ in range(12):zp.append(mul(zp[-1],z))
    for _ in range(4):wp.append(mul(wp[-1],w))
    val=(F(0),F(0)); az=val; aw=val
    for i,j,c in coeffs:
        val=add(val,scale(c,mul(zp[i],wp[j])))
        if i:az=add(az,scale(i*c,mul(zp[i-1],wp[j])))
        if j:aw=add(aw,scale(j*c,mul(zp[i],wp[j-1])))
    return val,az,aw

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate',nargs='?',type=Path,
                        default=Path(__file__).with_name('pure_4_12_certificate.json'))
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    raw=json.loads(args.certificate.read_text(encoding='utf-8'))
    coeffs=[(int(a['z']),int(a['w']),F(int(a['numerator']),int(a['denominator'])))
            for a in raw['coefficients']]
    assert len({(i,j) for i,j,c in coeffs})==len(coeffs),'Duplicate monomial'
    assert len(coeffs)==43
    kappa=F(int(raw['kappa'][0]),int(raw['kappa'][1]))
    tau=F(int(raw['tau'][0]),int(raw['tau'][1]))
    assert kappa>0 and tau>0
    assert {(i,j):c for i,j,c in coeffs}==reconstruct(kappa,tau),'Formula mismatch'
    assert all(c and 0<=i<=12 and 0<=j<=4 for i,j,c in coeffs)
    assert all(i+j<12 or (i,j)==(12,0) for i,j,c in coeffs)
    assert next(c for i,j,c in coeffs if (i,j)==(12,0))==1
    assert max(j for i,j,c in coeffs)==4
    den=int(raw['center_denominator'])
    r=F(int(raw['radius'][0]),int(raw['radius'][1]))
    assert den>0 and r>0
    centers=[(F(int(a),den),F(int(b),den)) for a,b in raw['centers']]
    assert len(centers)==144
    for i,(x,y) in enumerate(centers):
        for xx,yy in centers[i+1:]:
            assert abs(x-xx)>2*r or abs(y-yy)>2*r,'Overlapping boxes'
    assert sum(y==0 for x,y in centers)==12
    assert all(y==0 or abs(y)>r for x,y in centers),'Ambiguous real-axis box'

    signs=[]; qmax=F(0); beta_ratio_max=F(0); start=time.time()
    for index,(x,y) in enumerate(centers):
        (fr,fi),az,aw=evaluate(coeffs,x,y)
        a=az[0]+aw[0]; b=aw[1]-az[1]
        c=az[1]+aw[1]; d=az[0]-aw[0]
        det=a*d-b*c
        assert det,'Singular center Jacobian'
        inverse_norm=max(abs(d)+abs(b),abs(c)+abs(a))/abs(det)
        beta=max(abs(d*fr-b*fi),abs(-c*fr+a*fi))/abs(det)
        radius_bound=abs(x)+abs(y)+2*r
        powers=[F(1)]
        for _ in range(10):powers.append(powers[-1]*radius_bound)
        second_bound=sum((abs(co)*(i*(i-1)+i*j+j*(j-1))*powers[i+j-2]
                          for i,j,co in coeffs if i+j>=2),F(0))
        q=4*r*inverse_norm*second_bound
        assert q<1,('Contraction failed',index)
        assert beta+q*r<r,('Self-mapping failed',index)
        signs.append(1 if det>0 else -1)
        qmax=max(qmax,q); beta_ratio_max=max(beta_ratio_max,beta/r)
        if (index+1)%16==0:
            print('Verified',index+1,'of144 boxes; elapsed',round(time.time()-start,2),flush=True)
    assert signs.count(1)==78 and signs.count(-1)==66
    assert qmax<F(1,10**63)
    assert beta_ratio_max<F(1,10**60)
    scale_to_integers=lcm(*(c.denominator for i,j,c in coeffs))
    assert scale_to_integers==256*10**144
    report={
        'all_checks_passed':True,'arithmetic':'fractions.Fraction, exact throughout all assertions',
        'simple_real_slice_zeros':144,'positive_jacobian_zeros':78,'negative_jacobian_zeros':66,
        'real_axis_zeros':12,'disjoint_boxes':144,'pure_highest_part':'z^12',
        'analytic_degree':12,'conjugate_degree':4,
        'bezout_exhaustion':144,'polynomial_terms':len(coeffs),
        'independent_coefficient_reconstruction_passed':True,
        'kappa':[str(kappa.numerator),str(kappa.denominator)],
        'tau':[str(tau.numerator),str(tau.denominator)],
        'strict_exact_contraction_bound':'10^-63',
        'strict_exact_beta_over_radius_bound':'10^-60',
        'integer_scaling_denominator':str(scale_to_integers),
        'radius':[str(r.numerator),str(r.denominator)],
        'display_only_max_contraction_factor':float(qmax),
        'display_only_max_beta_over_radius':float(beta_ratio_max),
        'seconds':round(time.time()-start,2)}
    if args.report:args.report.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
