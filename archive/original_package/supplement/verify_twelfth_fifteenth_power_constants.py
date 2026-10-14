"""Exact rational certificates for C_12 and C_15. No sampled maximum is used."""
from pathlib import Path
from math import isqrt
import sympy as s,hashlib,json,time
if not __debug__:raise SystemExit('Assertions are required; do not run with -O.')
X,Y,R,U,V,a,b,B=s.symbols('X Y R U V a b B')

def newton(n,u,v):
    pp=[s.S(3),u,s.expand(u*u-2*v)]
    for j in range(3,n+1):pp.append(s.expand(u*pp[-1]-v*pp[-2]+pp[-3]))
    return pp[n]

def construct_q(k):
    t=a+s.I*b;tc=a-s.I*b
    P=newton(k//3,t,tc)-3
    real=s.S(0);imag=s.S(0)
    for (j,),c in s.Poly(P,b).terms():
        if j%2:imag+=c/s.I*B**((j-1)//2)
        else:real+=c*B**(j//2)
    ar=3-3*Y+Y*X;br=Y**2*(Y-X**2)
    rr=s.expand(real.subs({a:ar,B:br}));ii=s.expand(imag.subs({a:ar,B:br}))
    Q=s.cancel((rr*rr+br*ii*ii)/Y)
    assert s.Poly(Q,X,Y).domain==s.ZZ
    # Independent recurrence in the original trace verifies the power map.
    original=newton(k,U,V)
    third=newton(3,U,V);thirdc=newton(3,V,U)
    assert s.expand(original-newton(k//3,third,thirdc))==0
    return Q,original

def variations(seq,at):
    signs=[s.sign(p.eval(at)) for p in seq]
    signs=[z for z in signs if z]
    return sum(a!=b for a,b in zip(signs,signs[1:]))

def interval(poly,variables,bounds):
    def mul(a,b):
        vals=[x*y for x in a for y in b]
        return min(vals),max(vals)
    result=[s.S(0),s.S(0)]
    for powers,c in s.Poly(poly,*variables).terms():
        term=(c,c)
        for box,n in zip(bounds,powers):
            power=(s.S(1),s.S(1))
            for _ in range(n):power=mul(power,box)
            term=mul(term,power)
        result[0]+=term[0];result[1]+=term[1]
    return tuple(result)

def boundary_certificate(poly,right):
    degree=s.degree(poly,R)
    pending=[(s.S(0),right)];done=[]
    while pending:
        lo,hi=pending.pop()
        p=s.Poly(poly.subs(R,lo+(hi-lo)*R).expand(),R)
        cs=[sum(p.nth(j)*s.binomial(i,j)/s.binomial(degree,j) for j in range(i+1)) for i in range(degree+1)]
        if min(cs)>0:
            done.append(dict(interval=[str(lo),str(hi)],coefficients=list(map(str,cs)),minimum=str(min(cs))))
        else:
            assert hi-lo>right/2**16,'Boundary certificate not found'
            mid=(lo+hi)/2
            pending.extend([(lo,mid),(mid,hi)])
    return done

def verify(k):
    started=time.time()
    lower={12:s.Rational(23,2),15:s.Rational(57,4)}[k]
    right=6/lower;upper=right*right
    Q,P=construct_q(k)
    print(k,'constructed',flush=True)
    qx,qy=s.diff(Q,X),s.diff(Q,Y)
    resultant=s.Poly(s.resultant(qx,qy,X),Y)
    content,factors=s.factor_list(resultant)
    print(k,'factored resultant',flush=True)
    reconstructed=s.Poly(content,Y)
    root_rows=[];selected=[]
    for f,m in factors:
        reconstructed*=f**m
        seq=[s.Poly(p,Y) for p in s.sturm(f.as_expr(),Y)]
        v=[variations(seq,at) for at in [s.S(0),upper]]
        count=v[0]-v[1]
        # The Y factor has a root at zero, omitted by Sturm's endpoint convention.
        assert f.eval(upper)!=0
        root_rows.append(dict(polynomial=str(f.as_expr()),multiplicity=int(m),degree=f.degree(),
                              sturm_variations=v,roots_in_open_interval=count))
        if count:selected.append((f,count))
    assert reconstructed==resultant
    assert len(selected)==1 and selected[0][1]==1
    F=selected[0][0]
    assert F.degree()=={12:18,15:29}[k]
    print(k,'Sturm checks passed',flush=True)
    # The lexicographic basis is recomputed from the full gradient ideal.
    gb=s.groebner([qx,qy,F.as_expr()],X,Y,order='lex',domain=s.QQ)
    assert len(gb.polys)==2
    relation=gb.polys[0].as_expr()
    assert s.diff(relation,X)==1
    xp=s.expand(X-relation)
    assert not xp.has(X)
    assert gb==s.groebner([X-xp,F.as_expr()],X,Y,order='lex',domain=s.QQ)
    modulus=s.Poly(F,Y,domain=s.QQ)
    xpmod=s.Poly(xp,Y,domain=s.QQ)
    def reduced_substitution(poly):
        p=s.Poly(poly,X)
        answer=s.Poly(0,Y,domain=s.QQ)
        for c in p.all_coeffs():
            answer=(answer*xpmod+s.Poly(c,Y,domain=s.QQ)).rem(modulus)
        return answer
    assert reduced_substitution(qx).is_zero
    assert reduced_substitution(qy).is_zero
    print(k,'gradient ideal checked',flush=True)
    boundary=[]
    for sigma in [1,-1]:
        g=s.cancel((3-P.subs({U:sigma*R,V:sigma*R}))/R)
        assert s.expand(Q.subs({X:sigma*R,Y:R*R})-g*g)==0
        # On an admissible real trace, p_k=1+2 cos(k theta), hence g>=0.
        row=boundary_certificate(lower-g,right)
        boundary.append(dict(sign=sigma,polynomial=str(g),bernstein_pieces=row))
    print(k,'boundary checks passed',flush=True)
    lo,hi=s.S(0),upper
    slo=s.sign(F.eval(lo));assert slo and s.sign(F.eval(hi))==-slo
    for _ in range(150):
        mid=(lo+hi)/2
        if s.sign(F.eval(mid))==slo:lo=mid
        else:hi=mid
    xl,xh=interval(xp,[Y],[(lo,hi)])
    assert xl<xh<0 and xl*xl<lo<hi<upper<1
    vl,vh=interval(Q,[X,Y],[(xl,xh),(lo,hi)])
    assert vl>lower*lower
    den=10**10
    cl=s.Rational(isqrt(int(s.floor(vl*den**2))),den)
    ch=s.Rational(isqrt(int(s.floor(vh*den**2)))+1,den)
    assert cl*cl<vl<vh<ch*ch and ch-cl==s.Rational(1,den)
    print(k,'certified',str(s.N(cl,15)),str(s.N(ch,15)),'seconds',time.time()-started,flush=True)
    return dict(k=k,Q=str(Q),cutoff_R=str(right),boundary_threshold=str(lower),
                resultant_content=str(content),resultant_factors=root_rows,
                stationary_polynomial=str(F.as_expr()),stationary_degree=F.degree(),
                x_polynomial=str(xp),boundary_certificates=boundary,
                Y_interval=[str(lo),str(hi)],X_interval=[str(xl),str(xh)],
                C_interval=[str(cl),str(ch)],C_decimal=[str(s.N(cl,15)),str(s.N(ch,15))],
                equality='six traces: |s|^2=Y and Re(s^3)=XY',seconds=time.time()-started)

def main():
    start=time.time()
    rows=[verify(k) for k in [12,15]]
    report=dict(all_passed=True,arithmetic='exact rational polynomial arithmetic',
                orders=rows,verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                seconds=time.time()-start)
    target=Path(__file__).with_name('twelfth_fifteenth_power_constants_verification.json')
    target.write_text(json.dumps(report,indent=2)+'\n')
    print('ALL PASSED',target,flush=True)

if __name__=='__main__':main()
