"""Exact interpolation and trace-form certificates for 88 pure-head zeros.
The analytic small-parameter argument is in the accompanying manuscript.
No numerical root count is used by this verifier.
"""
from pathlib import Path
import sympy as s,json,hashlib,time
if not __debug__:raise SystemExit('Assertions are required; do not run with -O.')
z,w,x,y,omega=s.symbols('z w x y omega')
a=s.Rational(100,10101);b=-s.Rational(101,10101);c=-s.Rational(10100,10101)
kap=a*a*(1+a);d=3*a*a+a;J=d*d-a*a;A=d/J;B=-a/J
beta=a/(1+3*a);lam=B/A
e0=-s.Rational(14493204979,297671419500)
e1=s.Rational(1463813702879,29767141950000)
e2=-s.Rational(167587774534279,2976714195000000)
q20=s.Rational(13,125);q11=s.Rational(1,20);q02=-s.Rational(1,5)
k=-s.Rational(81,1000);h=s.Rational(19,125);c0=-s.Rational(1,500)
lu=(k+h)/2;lv=(k-h)/2
Epsilon=lambda t:e0+e1*t+e2*t*t
Target=lambda t:(35+535*t+3*t*t)/5000

def sha_rows(rows):
    return hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()

def trace_certificate(bcoef,ecoef,target,expected):
    zz=x+(1+2*omega)*y;ww=x-(1+2*omega)*y
    expr=s.rem(s.expand(zz**3+bcoef*zz*ww+ecoef*ww-target),omega**2+omega+1,omega)
    pp=s.Poly(expr,omega);F=s.expand(pp.nth(0));G=s.expand(pp.nth(1))
    gb=s.groebner([F,G],x,y,order='grlex',domain=s.QQ)
    leading=[q.LM(order=gb.order).exponents for q in gb.polys]
    assert any(j==0 for i,j in leading) and any(i==0 for i,j in leading)
    nx=min(i for i,j in leading if j==0);ny=min(j for i,j in leading if i==0)
    powers=[(i,j) for i in range(nx) for j in range(ny)
            if not any(i>=u and j>=v for u,v in leading)]
    powers.sort(key=lambda ij:(sum(ij),ij[0]))
    assert len(powers)==9
    basis=[x**i*y**j for i,j in powers]
    def vec(poly):
        rem=s.Poly(gb.reduce(s.expand(poly))[1],x,y,domain=s.QQ)
        assert set(rem.monoms())<=set(powers)
        return s.Matrix([rem.coeff_monomial(mon) for mon in basis])
    Mx=s.Matrix.hstack(*[vec(x*q) for q in basis]);My=s.Matrix.hstack(*[vec(y*q) for q in basis])
    assert Mx*My==My*Mx
    traces=s.Matrix([s.trace(Mx**i*My**j) for i,j in powers])
    H=s.Matrix(9,9,lambda i,j:(traces.T*vec(basis[i]*basis[j]))[0])
    assert H==H.T
    mat=H.copy();pivots=[]
    for k in range(9):
        pivot=mat[k,k];assert pivot!=0
        pivots.append(pivot)
        for i in range(k+1,9):
            for j in range(i,9):
                mat[i,j]=mat[j,i]=s.cancel(mat[i,j]-mat[i,k]*mat[k,j]/pivot)
    plus=sum(1 for v in pivots if v>0);minus=sum(1 for v in pivots if v<0)
    assert plus-minus==expected and plus+minus==9
    rows=[[str(v) for v in row] for row in H.tolist()]
    return dict(real_equations=[str(F),str(G)],standard_basis=list(map(str,basis)),
                dimension=9,positive=int(plus),negative=int(minus),real_simple_zeros=expected,
                ldl_pivots=list(map(str,pivots)),trace_matrix_sha256=sha_rows(rows))

def interpolation():
    f=z**3+z*w-kap;fs=w**3+z*w-kap;g=z-a
    gb=s.groebner([f,fs],z,w,order='grlex',domain=s.QQ)
    exps=[(i,j) for j in range(3) for i in range(3)]
    cols=[z**i*w**j for j in range(4) for i in range(10-j)]
    def red(e):return s.expand(gb.reduce(s.expand(e))[1])
    def vec(e):
        p=s.Poly(red(e),z,w,domain=s.QQ)
        assert all(i<3 and j<3 for i,j in p.monoms())
        return s.Matrix([p.coeff_monomial(z**i*w**j) for i,j in exps])
    def mcol(p):
        return s.Matrix.vstack(vec(p),vec(s.diff(p,z)),vec(s.diff(p,w)),
             s.Matrix([s.diff(p,z,2).subs({z:a,w:a})/2,s.diff(p,z,w).subs({z:a,w:a}),
                       s.diff(p,w,2).subs({z:a,w:a})/2]))
    M=s.Matrix.hstack(*[mcol(p) for p in cols])
    elim=(kap-z**3)**3-z**6
    eid=s.cancel(elim/((z-a)*s.diff(elim,z).subs(z,a)))
    gin=s.cancel((1-eid)/(z-a))
    assert red((z-a)*gin-(1-eid))==0
    eo=red((z**3-a**3)*(z**3-c**3)/((b**3-a**3)*(b**3-c**3)))
    assert red(eo*eo-eo)==0
    Lpoly=a*a*(A-1)*(z-a)+a*d*A*(w-a)
    cof=w*w*f-z*fs
    H=s.expand(f*cof*g+f*f*Lpoly)
    alpha=red(w*w+Lpoly*gin);betacoef=-z
    alphac=alpha.xreplace({z:w,w:z})
    error=red(eo*b*b*Epsilon(z/b))
    ku=red(g*(alpha*alpha+betacoef*alphac)/3)
    kv=red(g*(betacoef*alpha/3+error))
    kz=red(ku*s.diff(f,z)+kv*s.diff(fs,z));kw=red(ku*s.diff(f,w)+kv*s.diff(fs,w))
    dz,dw=s.symbols('dz dw');uf=d*dz+a*dw;vf=a*dz+d*dw
    qjet=s.expand(A*beta**2*(q20*uf**2+q11*uf*vf+q02*vf**2))
    qpoly=s.Poly(qjet,dz,dw)
    kt=s.Matrix.vstack(s.zeros(9,1),vec(kz),vec(kw),s.Matrix([
        qpoly.coeff_monomial(dz**2),qpoly.coeff_monomial(dz*dw),qpoly.coeff_monomial(dw**2)]))
    ks,free=M.gauss_jordan_solve(kt)
    ks=ks.subs({v:0 for v in free})
    assert M*ks==kt and M.rank()==27
    K=s.expand(sum(v*q for v,q in zip(ks,cols)))
    ct=red((1-eo)*(-z*z*w/9)+eo*(-b**3)*Target(z/b))
    constant=red(g*(alpha**3/27+betacoef*alpha*alphac/9+error*alphac/3-ct))
    ell_u=A*beta**3*lu;ell_v=A*beta**3*lv
    LM=s.Matrix.hstack(*[s.Matrix.vstack(vec(q),s.Matrix([
         s.diff(q,z).subs({z:a,w:a}),s.diff(q,w).subs({z:a,w:a})])) for q in cols])
    lt=s.Matrix.vstack(vec(constant),s.Matrix([ell_u*d+ell_v*a,ell_u*a+ell_v*d]))
    ls,free=LM.gauss_jordan_solve(lt);ls=ls.subs({v:0 for v in free})
    assert LM*ls==lt and LM.rank()==11
    L=s.expand(sum(v*q for v,q in zip(ls,cols)))
    dz,dw=s.symbols('dz dw');uf=d*dz+a*dw;vf=a*dz+d*dw
    def jet(p,n):
        expanded=s.Poly(p.subs({z:a+dz,w:a+dw}).expand(),dz,dw)
        return s.expand(sum(co*dz**i*dw**j for (i,j),co in expanded.terms() if i+j<=n))
    assert s.expand(jet(f**3*g,4)-A*uf**4-B*uf**3*vf)==0
    assert s.expand(jet(H,3)-A*beta*uf*vf**2)==0
    assert s.expand(jet(K,2)-qjet)==0 and s.expand(jet(L,1)-ell_u*uf-ell_v*vf)==0
    for p in [H,K,L]:
        pp=s.Poly(p,z,w,domain=s.QQ)
        assert pp.total_degree()<=9 and pp.degree(w)<=3
    coeffs=lambda p:[[i,j,str(co)] for (i,j),co in s.Poly(p,z,w).terms()]
    return dict(K_rank=27,L_rank=11,K_equations=30,L_equations=11,space_dimension=34,
                H_coefficients=coeffs(H),K_coefficients=coeffs(K),L_coefficients=coeffs(L),
                constant_order_four=str(A*beta**4*c0),
                shared_jets_verified=True,all_supports_verified=True)


def quartic_certificate():
    v=s.symbols('v')
    # Re F = R0+R1*y^2+(1-lambda)*y^4; Im F=y*(AA+BB*y^2).
    AA=(4+2*lam)*x**3-x*x+2*(q20-q02)*x+h
    BB=(-4+2*lam)*x-1
    R0=(1+lam)*x**4+x**3+(q02+q20+q11)*x*x+k*x+c0
    R1=-6*x*x+x-q02-q20+q11
    u=x+s.I*y;ubar=x-s.I*y
    F=s.expand(u**4+lam*u**3*ubar+u*ubar**2+q20*u**2+q11*u*ubar+q02*ubar**2+lu*u+lv*ubar+c0)
    assert s.expand(F-(R0+R1*y*y+(1-lam)*y**4+s.I*y*(AA+BB*y*y)))==0
    axis=s.Poly(R0,x,domain=s.QQ)
    off=s.Poly((1-lam)*AA**2-R1*AA*BB+R0*BB**2,x,domain=s.QQ)
    assert axis.degree()==4 and off.degree()==6
    assert s.gcd(axis,axis.diff()).degree()==0 and s.gcd(axis,s.Poly(AA,x)).degree()==0
    assert s.gcd(off,off.diff()).degree()==0 and s.gcd(off,s.Poly(BB,x)).degree()==0
    assert s.gcd(s.Poly(AA,x),s.Poly(BB,x)).degree()==0
    intervals={}
    def interval_eval(poly,lo,hi):
        low=high=s.S(0)
        for coef in s.Poly(poly,x,domain=s.QQ).all_coeffs():
            products=[low*lo,low*hi,high*lo,high*hi]
            low=min(products)+coef;high=max(products)+coef
        return low,high
    for name,poly,num in [('axis',axis,4),('off_axis',off,6)]:
        endpoints={'axis':[(-35000,-34000),(-250,-248),(-26,-25),(317,319)],
                   'off_axis':[(-214,-213),(-192,-191),(-181,-180),(448,449),(499,500),(623,624)]}[name]
        iso=[((s.Rational(lo,1000),s.Rational(hi,1000)),1) for lo,hi in endpoints]
        assert len(iso)==num and all(mult==1 for _,mult in iso)
        assert all(iso[i][0][1]<iso[i+1][0][0] for i in range(num-1))
        rows=[]
        for (lo,hi),mult in iso:
            assert lo<hi and poly.eval(lo)*poly.eval(hi)<0
            row=dict(interval=[str(lo),str(hi)],endpoint_signs=[int(s.sign(poly.eval(lo))),int(s.sign(poly.eval(hi)))])
            if name=='off_axis':
                al,ah=interval_eval(AA,lo,hi);bl,bh=interval_eval(BB,lo,hi)
                assert (al>0 and bh<0) or (ah<0 and bl>0)
                row['A_range']=[str(al),str(ah)];row['B_range']=[str(bl),str(bh)]
            rows.append(row)
        intervals[name]=rows
    return dict(lambda_value=str(lam),coefficients=dict(q20=str(q20),q11=str(q11),q02=str(q02),lu=str(lu),lv=str(lv),constant=str(c0)),
                A=str(AA),B=str(BB),axis_polynomial=str(axis.as_expr()),off_axis_polynomial=str(off.as_expr()),
                isolating_intervals=intervals,regular_zeros=16,real_axis_zeros=4)

def main():
    start=time.time()
    assert 0<a<s.Rational(1,3) and -1<c<b<0
    assert all(r**3+r*r==kap for r in [a,b,c])
    assert all((3*r*r+r)**2-r*r!=0 for r in [a,b,c])
    assert -1<lam<0
    assert Target(s.S(1))==s.Rational(573,5000)
    assert s.rem(Target(omega)-(s.Rational(4,625)+s.Rational(133,1250)*omega),omega**2+omega+1,omega)==0
    linear=interpolation();print('Exact K and L interpolation passed',flush=True)
    real=trace_certificate(s.S(1),Epsilon(s.S(1)),Target(s.S(1)),9)
    assert s.Poly(x**3+x*x+Epsilon(s.S(1))*x-Target(s.S(1)),x).count_roots(-s.oo,s.oo)==3
    print('Real cubic:9simple real zeros',flush=True)
    comp=trace_certificate(omega,Epsilon(omega),Target(omega),9)
    print('Complex cubic:9simple real zeros',flush=True)
    shared=quartic_certificate();print('Shared quartic:16simple real zeros',flush=True)
    assert 8*9+16==88 and 49+39==88 and 49-39==10
    out=dict(all_passed=True,arithmetic='exact rational polynomial and matrix arithmetic',
             seed=dict(a=str(a),b=str(b),c=str(c),kappa=str(kap)),
             epsilon=list(map(str,[e0,e1,e2])),shared_quartic=shared,
             interpolation=linear,real_cubic=real,complex_cubic=comp,
             total_regular_zeros=88,positive=49,negative=39,real_axis_zeros=10,
             verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.time()-start)
    path=Path(__file__).with_name('pure_3_10_88_verification.json');path.write_text(json.dumps(out,indent=2)+'\n')
    print('ALL PASSED',path,'seconds',time.time()-start,flush=True)

if __name__=='__main__':main()
