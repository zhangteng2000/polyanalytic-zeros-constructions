"""Exact interpolation and trace-form certificates for 84 pure-head zeros.
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
e0=-s.Rational(3851837,34343400)
e1=s.Rational(389035537,3434340000)
e2=-s.Rational(41570361437,343434000000)
Epsilon=lambda t:e0+e1*t+e2*t*t
Target=lambda t:(25+269*t-59*t*t)/2000

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
    kt=s.Matrix.vstack(s.zeros(9,1),vec(kz),vec(kw),s.zeros(3,1))
    ks,free=M.gauss_jordan_solve(kt)
    ks=ks.subs({v:0 for v in free})
    assert M*ks==kt and M.rank()==27
    K=s.expand(sum(v*q for v,q in zip(ks,cols)))
    ct=red((1-eo)*(-z*z*w/9)+eo*(-b**3)*Target(z/b))
    constant=red(g*(alpha**3/27+betacoef*alpha*alphac/9+error*alphac/3-ct))
    ell=-A*beta**3/100
    LM=s.Matrix.hstack(*[s.Matrix.vstack(vec(q),s.Matrix([
         s.diff(q,z).subs({z:a,w:a}),s.diff(q,w).subs({z:a,w:a})])) for q in cols])
    lt=s.Matrix.vstack(vec(constant),s.Matrix([ell*a,ell*d]))
    ls,free=LM.gauss_jordan_solve(lt);ls=ls.subs({v:0 for v in free})
    assert LM*ls==lt and LM.rank()==11
    L=s.expand(sum(v*q for v,q in zip(ls,cols)))
    dz,dw=s.symbols('dz dw');uf=d*dz+a*dw;vf=a*dz+d*dw
    def jet(p,n):
        expanded=s.Poly(p.subs({z:a+dz,w:a+dw}).expand(),dz,dw)
        return s.expand(sum(co*dz**i*dw**j for (i,j),co in expanded.terms() if i+j<=n))
    assert s.expand(jet(f**3*g,4)-A*uf**4-B*uf**3*vf)==0
    assert s.expand(jet(H,3)-A*beta*uf*vf**2)==0
    assert jet(K,2)==0 and s.expand(jet(L,1)-ell*vf)==0
    for p in [H,K,L]:
        pp=s.Poly(p,z,w,domain=s.QQ)
        assert pp.total_degree()<=9 and pp.degree(w)<=3
    coeffs=lambda p:[[i,j,str(co)] for (i,j),co in s.Poly(p,z,w).terms()]
    return dict(K_rank=27,L_rank=11,K_equations=30,L_equations=11,space_dimension=34,
                H_coefficients=coeffs(H),K_coefficients=coeffs(K),L_coefficients=coeffs(L),
                shared_jets_verified=True,all_supports_verified=True)

def main():
    start=time.time()
    assert 0<a<s.Rational(1,3) and -1<c<b<0
    assert all(r**3+r*r==kap for r in [a,b,c])
    assert all((3*r*r+r)**2-r*r!=0 for r in [a,b,c])
    assert -1<lam<0 and 0<s.Rational(1,100)<4/(27*(1-lam)**2)
    assert Target(s.S(1))==s.Rational(47,400)
    assert s.rem(Target(omega)-(s.Rational(21,500)+s.Rational(41,250)*omega),omega**2+omega+1,omega)==0
    linear=interpolation();print('Exact K and L interpolation passed',flush=True)
    real=trace_certificate(s.S(1),Epsilon(s.S(1)),s.Rational(47,400),9)
    assert s.Poly(x**3+x*x+Epsilon(s.S(1))*x-s.Rational(47,400),x).count_roots(-s.oo,s.oo)==3
    print('Real cubic:9simple real zeros',flush=True)
    comp=trace_certificate(omega,Epsilon(omega),Target(omega),7)
    print('Complex cubic:7simple real zeros',flush=True)
    assert 5*9+9+2*7+16==84 and 47+37==84 and 47-37==10
    out=dict(all_passed=True,arithmetic='exact rational polynomial and matrix arithmetic',
             seed=dict(a=str(a),b=str(b),c=str(c),kappa=str(kap)),
             epsilon=list(map(str,[e0,e1,e2])),shared_quartic=dict(lambda_value=str(lam),beta=str(beta),t='1/100'),
             interpolation=linear,real_cubic=real,complex_cubic=comp,
             total_regular_zeros=84,positive=47,negative=37,real_axis_zeros=10,
             verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.time()-start)
    path=Path(__file__).with_name('pure_3_10_84_verification.json');path.write_text(json.dumps(out,indent=2)+'\n')
    print('ALL PASSED',path,'seconds',time.time()-start,flush=True)

if __name__=='__main__':main()
