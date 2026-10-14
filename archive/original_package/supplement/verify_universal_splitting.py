"""Exact diagnostics for arbitrary-order splitting and the continuous sharp band.

Run with Python 3 and SymPy, without -O. No floating-point root count is used.
The all-parameter theorem is proved in the accompanying TeX/PDF; these finite
checks independently exercise its algebra and integer construction.
"""
from pathlib import Path
import sympy as s
import json, time

u,v,a,ac,b,bc=s.symbols('u v a ac b bc')
gens=(u,v,a,ac,b,bc)
z,w,x=s.symbols('z w x')
def local_terms(k):
    def p(e): return s.Poly(e,*gens,domain=s.QQ)
    U=p(u+a/s.Integer(k)); V=p(v+ac/s.Integer(k))
    beta=p(b/s.Integer(k)); alpha=p(bc/s.Integer(k))
    rec=[p(3),U,U*U+2*beta*V]
    for j in range(3,k+1):
        rec.append(U*rec[-1]+beta*V*rec[-2]-alpha*beta*beta*rec[-3])
    H=rec[k]+(-alpha)**k
    terms=[{} for _ in range(k+1)]
    for mon,c in H.terms():
        assert sum(mon)==k
        terms[k-mon[0]-mon[1]][mon]=c
    return terms

def p_coeffs(k,j):
    return [s.Rational(k,k-l)*s.binomial(k-l,l)*s.binomial(k-2*l,j-l)/s.Integer(k)**j
            for l in range(min(j,k-j)+1)]

def series_mul(A,B,L):
    return [sum(A[t]*B[h-t] for t in range(h+1)) for h in range(L+1)]

def error_coeffs(k,j,q):
    L=j-1
    inv=[s.S.One]+[s.S.Zero]*L
    for h in range(1,L+1): inv[h]=-sum(q[t]*inv[h-t] for t in range(1,h+1))
    power=q[:L+1]
    result=[s.S.Zero]*(L+1)
    for l,hc in enumerate(p_coeffs(k,j)):
        if l: power=series_mul(power,inv,L)
        result=[val+(-1)**l*hc*power[h] for h,val in enumerate(result)]
    assert result[0]==0
    return result

def verify(k,r):
    start=time.time(); N=2*r+1
    def p(e): return s.Poly(e,z,w,domain=s.QQ_I)
    lam=1+2*s.I
    Q=p(1)
    for h in range(1,r+1): Q*=p(z*w-h*h)
    q=[Q.coeff_monomial(z**(r-h)*w**(r-h)) for h in range(r+1)]
    f=p(z**N)+Q; fs=p(w**N)+Q
    polys=[f,fs,p(lam*w**(r+1)),p(s.conjugate(lam)*z**(r+1)),p(-lam*z**r),p(-s.conjugate(lam)*w**r)]
    powers=[[p(1)] for _ in polys]
    for j in range(6):
        for t in range(k): powers[j].append(powers[j][-1]*polys[j])
    raw=[]
    for terms in local_terms(k):
        out=p(0)
        for mon,c in terms.items():
            term=p(c)
            for l,h in enumerate(mon):
                if h: term*=powers[l][h]
            out+=term
        raw.append(out)
    assert raw[0]==f**k
    result=[]
    for j in range(1,k):
        hc=p_coeffs(k,j)
        assert sum((-1)**l*val for l,val in enumerate(hc))==0
        coeff=error_coeffs(k,j,q)
        C=p(0)
        for h in range(1,j):
            C+=p(lam**j*coeff[h]*z**(r-h)*w**((j-1)*r+j-1-h))
        T=raw[j]-C*powers[0][k-j-1]*fs
        assert T.degree(w)<=k*r,(k,r,j,T.degree(w),k*r)
        assert T.total_degree() < k*N
        result.append([j,int(raw[j].degree(w)),int(T.degree(w)),int(T.total_degree())])
    T=raw[k].reorder(w,z).rem(fs.reorder(w,z)).reorder(z,w)
    assert T.degree(w)<N and T.degree(w)<=k*r
    assert T.total_degree()<k*N
    assert (raw[k]-T).reorder(w,z).rem(fs.reorder(w,z)).is_zero
    result.append([k,int(raw[k].degree(w)),int(T.degree(w)),int(T.total_degree())])
    report={'k':k,'r':r,'degrees':result,'seconds':round(time.time()-start,3)}
    print(json.dumps(report),flush=True)
    return report


def cyclotomic_model_check(k):
    # Independent ring calculation: s^k=1, t^k=-1, alpha and beta formal.
    def mon(es,et,ea,eb,c=1):
        return {(es%k,et%k,ea,eb):c*(-1 if (et//k)%2 else 1)}
    def add(*values):
        result={}
        for value in values:
            for m,c in value.items(): result[m]=result.get(m,0)+c
        return {m:c for m,c in result.items() if c}
    def scale(value,c): return {m:c*v for m,v in value.items() if c*v}
    def mul(A,B):
        result={}
        for ma,ca in A.items():
            for mb,cb in B.items():
                m=tuple(a+b for a,b in zip(ma,mb))
                value=mon(*m,ca*cb)
                for n,c in value.items(): result[n]=result.get(n,0)+c
        return {m:c for m,c in result.items() if c}
    U=add(mon(-1,0,0,1),mon(0,-1,0,1),mon(1,1,1,0,-1))
    V=add(mon(1,0,1,0),mon(0,1,1,0),mon(-1,-1,0,1,-1))
    betaV=mul(mon(0,0,0,1),V)
    alpha_beta2=mon(0,0,1,2)
    rec=[mon(0,0,0,0,3),U,add(mul(U,U),scale(betaV,2))]
    for h in range(3,k+1):
        rec.append(add(mul(U,rec[-1]),mul(betaV,rec[-2]),scale(mul(alpha_beta2,rec[-3]),-1)))
    assert not add(rec[k],mon(0,0,k,0,(-1)**k))
    for es,et in [(1,0),(0,1)]:
        assert not add(mon(3*es,3*et,1,0),scale(mul(V,mon(2*es,2*et,0,0)),-1),
                       scale(mul(U,mon(es,et,0,0)),-1),mon(0,0,0,1))

def arithmetic_parameters(h):
    from math import isqrt
    S=isqrt(4*h)
    if S*S<4*h: S+=1
    K=S//2; delta=S-2*K; T=K*(K+delta)-h
    j=(isqrt(delta*delta+4*T)-delta)//2
    a=T-j*(j+delta); k=K-j; l=K+delta+j
    assert 0<=T<=K-1
    assert j*(j+delta)<=T<(j+1)*(j+1+delta)
    assert 0<=a<=2*j+delta and h+a==k*l and k<=l
    assert k-2*a>=K-5*j-2*delta
    return k,l,a

def run():
    if not __debug__: raise RuntimeError('Assertions must be enabled; do not use -O.')
    for k in range(2,13): cyclotomic_model_check(k)
    print('Independent cyclotomic model and recovery identities: 2 <= k <= 12 passed.',flush=True)
    p_checks=0
    for k in range(2,41):
        for j in range(1,k):
            assert sum((-1)**l*h for l,h in enumerate(p_coeffs(k,j)))==0
            p_checks+=1
    weighted_checks=0
    cs=s.symbols('c1:8')
    for k in range(2,9):
        for j in range(2,k):
            vals=error_coeffs(k,j,[s.S.One,*cs])
            for h in range(1,j):
                pp=s.Poly(s.expand(vals[h]),*cs,domain=s.QQ)
                for mon,c in pp.terms():
                    if c: assert sum((index+1)*power for index,power in enumerate(mon))==h
                weighted_checks+=1
    print('Intermediate cancellation and weighted-error identities passed.',flush=True)
    degree_records=[]
    for k in range(2,9): degree_records.append(verify(k,k-1))
    for k in range(2,7): degree_records.append(verify(k,k))
    arithmetic_count=250000
    for h in range(1,arithmetic_count+1): arithmetic_parameters(h)
    # Pure-degree padding + transfer + circle degree identities, symbolically.
    kk,ll,aa,dd=s.symbols('k l a d')
    hh=kk*ll-aa; ss=kk-dd-2*aa
    assert s.expand(kk*(ll-1)+aa+ss-(hh-dd))==0
    assert s.expand(kk*(2*ll-1)+ss-(2*hh-dd))==0
    fourth_power_pairs=0
    for t in range(1,41):
        m=t**4
        for d in range(1,max(0,t*t-5*t-4)+1):
            k,l,a=arithmetic_parameters(m+d)
            assert 2<=k<=l and d<=k-2*a
            assert k*(l-1)+a+(k-d-2*a)==m
            assert k*(2*l-1)+(k-d-2*a)==2*m+d
            fourth_power_pairs+=1
    assert arithmetic_parameters(101)==(8,13,3)  # The generic choice is sufficient, not optimal.
    assert 6*17==101+1 and 4<=6-2*1
    assert 97**2+198**2==48613
    assert 20**2+45**2==2425 and 22**2+47**2==2693
    assert 30**2+66**2==5256
    for d in range(1,47):
        k,l,a=arithmetic_parameters(10000+d)
        assert d<=k-2*a
    ratio=s.simplify((1+2*s.I)/(1-2*s.I))
    assert s.simplify(ratio+1/ratio)==-s.Rational(6,5)
    report={
        'status':'passed',
        'arithmetic':'exact integers, rationals, and Gaussian rationals; no floating-point root count',
        'cyclotomic_model_orders':list(range(2,13)),
        'intermediate_coefficient_identities':p_checks,
        'weighted_error_coefficient_checks':weighted_checks,
        'degree_cases':len(degree_records),
        'degree_specialization':'epsilon=1, solely for polynomial degree identities',
        'integer_construction_range':[1,arithmetic_count],
        'fourth_power_band_pairs':fourth_power_pairs,
        'm10000_band_checked':[1,46],
        'phase_ratio_trace':'-6/5, not a rational algebraic integer',
        'scope':'Finite diagnostics supplement the analytic infinite-family and band proofs. No fixed numerical small parameters are certified. The sufficient band boundary is not claimed optimal.',
        'cases':degree_records,
    }
    target=Path(__file__).with_name('universal_splitting_verification.json')
    target.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({key:val for key,val in report.items() if key!='cases'},indent=2),flush=True)

if __name__=='__main__': run()
