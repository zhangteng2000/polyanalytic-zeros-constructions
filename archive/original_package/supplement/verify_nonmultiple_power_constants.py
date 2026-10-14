"""Exact certificate for all nonmultiple-of-three power constants.
The proof of the infinite tail bound is supplied in the manuscript."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial
import time,json,hashlib
import sympy as s

H,U,V,X=s.symbols('H U V X')
K=36;CUT=20

def power_bernstein(poly,variables,intervals):
    q=s.Poly(s.expand(poly.subs({v:a+(b-a)*v for v,(a,b) in zip(variables,intervals)},simultaneous=True)),*variables)
    assert not q.as_expr().has(s.Float)
    assert all(c.is_Rational for _,c in q.terms())
    dims=[d+1 for d in q.degree_list()]
    coeff={ij:F(int(c.p),int(c.q)) for ij,c in q.terms()}
    for axis,dim in enumerate(dims):
        groups={}
        for idx,c in coeff.items():
            key=idx[:axis]+idx[axis+1:]
            groups.setdefault(key,{})[idx[axis]]=c
        out={}
        for key,row in groups.items():
            for j in range(dim):
                value=sum((c*F(comb(j,l),comb(dim-1,l)) for l,c in row.items() if l<=j),F(0))
                if value:out[key[:axis]+(j,)+key[axis:]]=value
        coeff=out
    expected=1
    for d in dims:expected*=d
    minimum=min(coeff.values()) if len(coeff)==expected else min(F(0),*coeff.values())
    digest=hashlib.sha256('\n'.join(str(coeff.get(idx,F(0))) for idx in __import__('itertools').product(*(range(d) for d in dims))).encode()).hexdigest()
    return dict(degrees=[d-1 for d in dims],minimum=str(minimum),min_float=float(minimum),count=expected,sha256=digest)

def truncated(residue):
    shift=0 if residue==1 else 2
    terms={}
    for d in range(1,CUT+1):
        for j in range(d+1):
            i=d-j
            if (i+2*j-residue)%3:continue
            exponent=(i-j-1+shift)//3
            ch=s.Poly(s.prod(1+s.Rational(3*l-i-2*j,K)*H for l in range(1,d)),H)
            for (a,),c in ch.terms():terms[exponent,d-1,a]=(-1)**j*c/s.factorial(i)/s.factorial(j)
    lo=min(t[0] for t in terms);hi=max(t[0] for t in terms)
    p=s.Poly.from_dict({(a-lo,b,c):v for (a,b,c),v in terms.items()},V,U,H)
    rev=s.Poly.from_dict({(hi-a,b,c):v for (a,b,c),v in terms.items()},V,U,H)
    product=p*rev;center=hi-lo
    pc=product.as_dict();maximum=max(t[0] for t in pc)-center
    phase=1 if residue==1 else -1
    dd={}
    for (vp,up,hp),c in pc.items():
        j=vp-center
        if j<=0:continue
        assert pc.get((center-j,up,hp),0)==c
        assert up>=1
        cp=s.Poly(s.diff(s.chebyshevt(j,X),X),X)
        for (xp,),a in cp.terms():dd[hp,up-1,xp]=dd.get((hp,up-1,xp),0)+phase*s.Rational(2,3)*a*c
    D=s.Poly.from_dict(dd,H,U,X).as_expr()
    radial={}
    for (vp,up,hp),c in terms.items():radial[hp,up]=radial.get((hp,up),0)+phase*s.Integer(phase)**vp*c
    G=s.Poly.from_dict(radial,H,U).as_expr()
    print('constructed',residue,'A terms',len(terms),'D terms',len(dd),flush=True)
    angular=power_bernstein(D,[H,U,X],[(s.S(0),s.S(1)),(s.S(0),s.S(1)),(-s.S(1),s.S(1))])
    concavity=power_bernstein(-s.diff(G,U,2)/9,[H,U],[(s.S(0),s.S(1)),(s.S(0),s.S(1))])
    end=power_bernstein(-s.diff(G,U).subs(U,1)/3,[H],[(s.S(0),s.S(1))])
    out=dict(residue=residue,angular=angular,concavity=concavity,negative_endpoint_derivative=end)
    print(out,flush=True);return out

def moment_bounds():
    d0=CUT+1
    def term(d):return F(2**d,factorial(d))*F(K+2*d,K)**(d-1)
    # e<3 gives this upper bound for successive absolute-degree terms.
    ratio=F(6,d0+1)+F(12,K)
    tails=[];totals=[]
    for j in range(3):
        rr=ratio*F(d0+1,d0)**j
        assert rr<1
        tail=F(d0**j)*term(d0)/(1-rr)
        tails.append(tail)
        totals.append(sum((d**j*term(d) for d in range(1,d0)),F(0))+tail)
    # Differentiating Chebyshev factors costs at most ((d+e)/3)^2;
    # division by r contributes 1/3 after bounding degree>=1 at r=3.
    angular=F(2,27)*(totals[2]*tails[0]+2*totals[1]*tails[1]+totals[0]*tails[2])
    concavity=tails[2]/9
    first=tails[1]/3
    return dict(tails=[str(v) for v in tails],totals=[str(v) for v in totals],
                angular_error=str(angular),concavity_error=str(concavity),endpoint_error=str(first),
                error_floats=[float(angular),float(concavity),float(first)])


def finite_order(k):
    z,w,R,t,v,x=s.symbols('z w R t v x')
    powers=[s.S(3),z,z*z-2*w]
    for j in range(3,k+1):powers.append(s.Poly(z*powers[-1]-w*powers[-2]+powers[-3],z,w).as_expr())
    sigma=1 if k%3==1 else -1
    shift=0 if sigma==1 else 2
    expression=s.expand(powers[k].subs({z:R*t,w:R/t})/(R*t)*t**shift)
    angular_terms={}
    for term in s.Add.make_args(expression):
        exponent=int(term.as_powers_dict().get(t,0));assert exponent%3==0
        j=exponent//3
        angular_terms[j]=angular_terms.get(j,0)+term/t**exponent
    A=s.expand(sum(c*v**j for j,c in angular_terms.items()))
    radial=s.cancel(powers[k].subs({z:sigma*R,w:sigma*R})/R)
    assert s.expand(sigma*A.subs(v,sigma)-radial)==0
    L=s.expand(A*A.subs(v,1/v));Q=L.coeff(v,0)
    for j in range(1,max(angular_terms)-min(angular_terms)+1):
        assert s.expand(L.coeff(v,j)-L.coeff(v,-j))==0
        Q+=2*L.coeff(v,j)*s.chebyshevt(j,x)
    assert s.expand(Q.subs(x,sigma)-radial**2)==0
    D=s.cancel(sigma*s.diff(Q,x)/R)
    right=s.Rational(3,k)
    angular=power_bernstein(D,[R,x],[(s.S(0),right),(-s.S(1),s.S(1))])
    concavity=power_bernstein(-s.diff(radial,R,2),[R],[(s.S(0),right)])
    assert F(angular['minimum'])>0 and F(concavity['minimum'])>0
    derivative=s.diff(radial,R)
    assert radial.subs(R,0)==k
    assert derivative.subs(R,0)>0>derivative.subs(R,right)
    assert abs(radial.subs(R,right))<=k
    lo,hi=s.S(0),right
    for _ in range(110):
        mid=(lo+hi)/2
        if derivative.subs(R,mid)>0:lo=mid
        else:hi=mid
    lower=upper=s.S(0)
    for (j,),a in s.Poly(radial,R).terms():
        aa,bb=lo**j,hi**j
        if a<0:aa,bb=bb,aa
        lower+=a*aa;upper+=a*bb
    assert lower>k
    den=10**10
    low=s.floor(lower*den)/den;high=s.ceiling(upper*den)/den
    assert high-low==s.Rational(1,den)
    out=dict(order=k,radial_polynomial=str(radial),angular=angular,concavity=concavity,
             root_interval=[str(lo),str(hi)],constant_interval=[str(low),str(high)],
             approximate_interval=[str(s.N(low,13)),str(s.N(high,13))])
    print('finite order',k,out['approximate_interval'],flush=True)
    return out

def independently_check_coefficient_formula():
    # The first computation is integer Newton recurrence. The second is the
    # product formula used for the uniform Taylor certificate.
    polys=[{(0,0):3},{(1,0):1},{(2,0):1,(0,1):-2}]
    for k in range(3,81):
        pp={}
        for prev,di,dj,sgn in [(polys[-1],1,0,1),(polys[-2],0,1,-1),(polys[-3],0,0,1)]:
            for (i,j),a in prev.items():pp[i+di,j+dj]=pp.get((i+di,j+dj),0)+sgn*a
        polys.append({ij:a for ij,a in pp.items() if a})
    checked=0
    for k in range(4,81):
        if k%3==0:continue
        for d in range(1,CUT+1):
            for j in range(d+1):
                i=d-j
                if (i+2*j-k)%3:continue
                product=F((-1)**j,factorial(i)*factorial(j))
                for l in range(1,d):product*=1+F(3*l-i-2*j,k)
                direct=F(polys[k].get((i,j),0),k)*F(3,k)**(d-1)
                assert product==direct,(k,i,j)
                checked+=1
    return checked

def main():
    start=time.time()
    independent=independently_check_coefficient_formula()
    finite=[finite_order(k) for k in range(10,36) if k%3]
    tail=moment_bounds()
    assert F(tail['angular_error'])<F(1,4000)
    assert F(tail['concavity_error'])<F(1,30000)
    assert F(tail['endpoint_error'])<F(1,200000)
    uniform=[truncated(r) for r in [1,2]]
    for row in uniform:
        assert F(row['angular']['minimum'])>F(1,9)
        assert F(row['concavity']['minimum'])>F(3,40)
        assert F(row['negative_endpoint_derivative']['minimum'])>F(1,6)
    expected=[['707efc139a4d3d9a3664cd18e7de2bf135f12a7071f2d70b7454280f3830b282',
               '1b8ddd64ca0d4a50a9389f81ef69af2cf2475aa31d565366aa418cf175deb621',
               '2bcf3e184d6e270393270b057a01c7ac78342e1aa898873f5e407920874ebd90'],
              ['4fe6629eba2944431c1a4e9d775be3cbefb0426a49eb8f8a117bf6e977dd9f95',
               '9f0f3e6412f726fdca3c2d33f04e7c85ce0d2c8e6dee39ba16be82c2457643fd',
               '76b67982951118e9de652048f253e4acb777ad74ba32b82ac6876654568ca255']]
    for row,hashes in zip(uniform,expected):
        assert [row[key]['sha256'] for key in ['angular','concavity','negative_endpoint_derivative']]==hashes
    out=dict(all_passed=True,arithmetic='rational polynomial arithmetic only',
             scope='all k >= 4 with k not divisible by 3; low orders 4,5,7,8 use the existing proofs',
             finite_orders=finite,k_lower=K,cutoff=CUT,tail=tail,uniform=uniform,
             independent_recurrence_checks=independent,seconds=time.time()-start,
             verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    target=Path(__file__).with_name('nonmultiple_power_constants_verification.json')
    target.write_text(json.dumps(out,indent=2)+'\n')
    print('ALL PASSED',target,flush=True)

if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions are required; do not run with -O.')
    main()
