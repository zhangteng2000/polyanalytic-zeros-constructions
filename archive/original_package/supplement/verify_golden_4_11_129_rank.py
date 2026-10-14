"""Exact finite-field witness for the 129-zero paired-seed construction.

This establishes nonvanishing of specified algebraic minors, not a floating
zero count. The analytic argument using these minors is documented separately.
All arithmetic below is integer arithmetic modulo a prime.
"""
from pathlib import Path
import json
import hashlib
from math import isqrt

P=0
N=6
def inv(x):return pow(x%P,-1,P)
def add(*polys):
    out={}
    for p in polys:
        for k,v in p.items():out[k]=(out.get(k,0)+v)%P
    return {k:v for k,v in out.items() if v}
def scale(p,c):return {k:v*c%P for k,v in p.items() if v*c%P}
def mul(p,q):
    out={}
    for (i,j),a in p.items():
        for (k,l),b in q.items():
            if i+j+k+l<=N:out[i+k,j+l]=(out.get((i+k,j+l),0)+a*b)%P
    return {k:v for k,v in out.items() if v}
def power(p,k):
    out={(0,0):1}
    for _ in range(k):out=mul(out,p)
    return out
def compose(p,z,w):
    if not p:return {}
    zs=[power(z,i) for i in range(max(i for i,j in p)+1)]
    ws=[power(w,j) for j in range(max(j for i,j in p)+1)]
    return add(*(scale(mul(zs[i],ws[j]),a) for (i,j),a in p.items()))
def star(p):return {(j,i):v for (i,j),v in p.items()}
def evaluate(p,z,w):return sum(a*pow(z,i,P)*pow(w,j,P) for (i,j),a in p.items())%P

def inverse_coordinates(p,z0,w0):
    h=compose(p,{(0,0):z0,(1,0):1},{(0,0):w0,(0,1):1})
    hs=compose(star(p),{(0,0):z0,(1,0):1},{(0,0):w0,(0,1):1})
    assert h.get((0,0),0)==hs.get((0,0),0)==0
    a,b=h.get((1,0),0),h.get((0,1),0)
    cc,d=hs.get((1,0),0),hs.get((0,1),0)
    det=(a*d-b*cc)%P;assert det
    Ji=[d*inv(det)%P,-b*inv(det)%P,-cc*inv(det)%P,a*inv(det)%P]
    dz={(1,0):Ji[0],(0,1):Ji[1]};dw={(1,0):Ji[2],(0,1):Ji[3]}
    hn={k:v for k,v in h.items() if sum(k)>=2};hsn={k:v for k,v in hs.items() if sum(k)>=2}
    for degree in range(2,N+1):
        r1=compose(hn,dz,dw);r2=compose(hsn,dz,dw)
        for i in range(degree+1):
            k=(i,degree-i);v1=r1.get(k,0);v2=r2.get(k,0)
            dz[k]=(-Ji[0]*v1-Ji[1]*v2)%P;dw[k]=(-Ji[2]*v1-Ji[3]*v2)%P
    z=add(dz,{(0,0):z0});w=add(dw,{(0,0):w0})
    assert compose(p,z,w)=={(1,0):1}
    assert compose(star(p),z,w)=={(0,1):1}
    return z,w,det

def primes(start,stop):
    for p in range(start,stop):
        if all(p%d for d in range(2,isqrt(p)+1)):yield p

def seed(prime):
    global P
    P=prime
    sq={a*a%P:a for a in range(P)}
    if P in [2,5] or -1%P not in sq or 5 not in sq:return None
    ii=sq[-1%P];c=inv(10)
    for root5 in [sq[5],-sq[5]%P]:
        ell=(root5-1)*inv(2)%P
        ds=(ell*ell-4*ell*c)%P
        gy2=(c*inv(1-ell)-inv(4))%P
        qd=((4-2*ell)**2-16*c*(3+ell))%P
        if any(v not in sq or not v for v in [ds,gy2,qd]):continue
        ax=[(-ell+s*sq[ds])*inv(2)%P for s in [-1,1]]
        off=[((-inv(2)+s*ii*sq[gy2])%P,(-inv(2)-s*ii*sq[gy2])%P,'off') for s in [-1,1]]
        nodes=[(a,a,'axis') for a in ax]+off+[(0,0,'fonly')]
        valid=True
        for s in [-1,1]:
            x=(-(4-2*ell)+s*sq[qd])*inv(8)%P
            yy=(3*x*x+2*x+c*(2+ell))%P
            if yy not in sq or not yy:valid=False;break
            nodes += [((x+t*ii*sq[yy])%P,(x-t*ii*sq[yy])%P,'fonly') for t in [-1,1]]
        if not valid or len({(z,w) for z,w,k in nodes})!=9:continue
        g={(2,0):1,(1,1):ell,(1,0):1,(0,0):c}
        f={(3,0):1,(2,0):1,(1,1):(ell-1)%P,(1,0):c*(1+ell)%P,(0,1):-c%P}
        if any(evaluate(f,z,w) or evaluate(star(f),z,w) for z,w,k in nodes):continue
        try:
            for z,w,k in nodes:inverse_coordinates(g if k!='fonly' else f,z,w)
        except (AssertionError,ValueError):continue
        return ell,c,ii,nodes,f,g
    return None

def reduce_matrix(A,b=None):
    A=[row[:] + ([] if b is None else [b[i]]) for i,row in enumerate(A)]
    n=len(A[0])-(b is not None);rank=0;pivots=[];pivots_raw=[]
    for col in range(n):
        ix=next((i for i in range(rank,len(A)) if A[i][col]%P),None)
        if ix is None:continue
        A[rank],A[ix]=A[ix],A[rank]
        v=A[rank][col]%P;pivots_raw.append(v);iv=inv(v)
        A[rank]=[v*iv%P for v in A[rank]]
        for i in range(len(A)):
            if i==rank:continue
            v=A[i][col]%P
            if v:A[i]=[(x-v*y)%P for x,y in zip(A[i],A[rank])]
        pivots.append(col);rank+=1
        if rank==len(A):break
    if b is None:return dict(rank=rank,pivots=pivots,pivots_raw=pivots_raw)
    assert all(not row[n] for row in A[rank:])
    x=[0]*n
    for row,col in enumerate(pivots):x[col]=A[row][n]
    return dict(rank=rank,pivots=pivots,pivots_raw=pivots_raw,solution=x)

def determinant(A):
    A=[row[:] for row in A];out=1
    assert all(len(row)==len(A) for row in A)
    for j in range(len(A)):
        k=next((k for k in range(j,len(A)) if A[k][j]%P),None)
        if k is None:return 0
        if k!=j:A[j],A[k]=A[k],A[j];out=-out%P
        v=A[j][j]%P;out=out*v%P;iv=inv(v)
        for k in range(j+1,len(A)):
            v=A[k][j]*iv%P
            for l in range(j+1,len(A)):A[k][l]=(A[k][l]-v*A[j][l])%P
            A[k][j]=0
    return out

def main():
    global P,N
    selected=None
    for prime in primes(101,10000):
        selected=seed(prime)
        if selected:break
    assert selected is not None
    ell,c,ii,nodes,f,g=selected
    support=[(i,j) for i in range(12) for j in range(5)]
    HA=[];hb=[];basehb=[];KA=[];CA=[];data=[]
    for idx,(z0,w0,kind) in enumerate(nodes):
        N=6
        z,w,det=inverse_coordinates(g if kind!='fonly' else f,z0,w0)
        mons=[mul(power(z,i),power(w,j)) for i,j in support]
        q={}
        if kind!='fonly':
            cube=power(compose(f,z,w),3)
            q={(a,b-1):v for (a,b),v in cube.items() if b>=1}
        hkeys=[(0,0),(1,0),(0,1)]
        if kind!='fonly':hkeys += [(2,0),(1,1),(0,2)]
        if kind=='axis':hkeys += [(d,0) for d in [3,4,5]]+[(0,d) for d in [3,4,5]]
        correction=0
        if kind=='off':
            LL=(z0-w0)%P;BB=ell*z0%P;ss=z0*w0%P
            delta=-pow(LL,3,P)*(3+LL*LL)*inv(4*ell*ell*ss*ss)%P
            correction=pow(BB,3,P)*delta%P
        for ij in hkeys:
            HA.append([p.get(ij,0) for p in mons])
            hb.append((q.get(ij,0)+(correction if ij==(2,0) else 0))%P)
            basehb.append(q.get(ij,0))
        kkeys=[(0,0),(1,0),(0,1)]
        if kind=='axis':kkeys += [(2,0),(0,2),(3,0),(0,3),(1,1)]
        for ij in kkeys:KA.append([p.get(ij,0) for p in mons])
        ckeys=[(0,0)]
        if kind=='axis':ckeys += [(1,0),(0,1)]
        for ij in ckeys:CA.append([p.get(ij,0) for p in mons])
        data.append(dict(z=z0,w=w0,kind=kind,det=det,mons=mons,q=q))
    hc=[int(ij==(11,4)) for ij in support]
    hs=reduce_matrix(HA,hb)
    assert hs['rank']==51
    initial=hs['solution']
    prescribed=initial[54]
    hs=reduce_matrix(HA,[(v-prescribed*row[58])%P for row,v in zip(HA,hb)])
    H=hs['solution'];H[58]=prescribed
    assert all((sum(a*x for a,x in zip(row,H))-v)%P==0 for row,v in zip(HA,hb))
    corner=sum(a*x for a,x in zip(hc,H))%P
    print('H corner',corner,'corner augmented rank',reduce_matrix(HA+[hc])['rank'],flush=True)
    assert corner==0
    assert H[58]==H[54]
    assert H[54]!=0
    nonzero_inner=[]
    for d in data[:2]:
        coeff=(d['q'].get((2,1),0)-sum(p.get((2,1),0)*x for p,x in zip(d['mons'],H)))%P
        assert coeff;nonzero_inner.append(coeff)
    kr=reduce_matrix(KA);cr=reduce_matrix(CA)
    assert kr['rank']==37 and cr['rank']==13
    minors={name:determinant([[row[j] for j in result['pivots']] for row in matrix]) for name,matrix,result in [('H',HA,hs),('K',KA,kr),('C',CA,cr)]}
    assert all(minors.values())
    assert (P,ell,c)==(461,21,415)
    assert minors=={'H':79,'K':260,'C':350}
    assert nonzero_inner==[247,251] and H[54]==H[58]==349
    print('modified axis coefficients',nonzero_inner,'linear infinity coefficient',H[54],flush=True)
    # The same selected K minor permits a prescribed corner coefficient 1.
    assert 59 not in kr['pivots']
    kc=[int(ij==(11,4)) for ij in support]
    assert reduce_matrix(KA+[kc])['rank']==38
    z0,w0,_=nodes[2];LL=(z0-w0)%P;BB=ell*z0%P
    DD=(ell*ell*z0*w0+LL*LL)%P
    infinity_linear=-8*(1+LL*LL)*inv(pow(ell,5,P)*DD*DD*pow(1-LL*LL,2,P))%P
    assert H[54]==infinity_linear
    J=reduce_matrix(HA,basehb)['solution']
    ss=z0*w0%P
    original_corner=-(1-ss)*inv(pow(ell,5,P)*DD*DD*pow(ss,3,P))%P
    assert J[59]==original_corner==47
    assert all((sum(a*x for a,x in zip(row,J))-v)%P==0 for row,v in zip(HA,basehb))
    # An explicit homogeneous lift changes the mixed cubic coefficient at
    # every f-only node while leaving the barred-square coefficient fixed.
    N=20
    zp={(1,0):1};wp={(0,1):1}
    qq={(1,0):(ell-1)%P,(0,0):-c%P}
    tt={(2,0):1,(1,0):ell,(0,0):ell*c%P}
    raw=scale(mul(mul(zp,qq),{(1,0):(2+ell)%P,(0,0):1}),ell*ell)
    ap=[raw.get((i,0),0) for i in range(4)]
    for degree in [3,2]:
        coeff=ap[degree];ap[degree]=0
        ap[degree-1]=(ap[degree-1]-coeff*ell)%P
        ap[degree-2]=(ap[degree-2]-coeff*ell*c)%P
    aa={(i,0):ap[i] for i in [0,1] if ap[i]}
    GG=add(mul(power(wp,2),f),scale(mul(qq,star(f)),-1))
    H0=mul(mul(g,f),add(mul(tt,GG),mul(aa,f)))
    assert max(i for i,j in H0)<=11 and max(j for i,j in H0)<=4 and max(i+j for i,j in H0)<=12
    N=6
    H0vec=[H0.get(ij,0) for ij in support]
    assert all(sum(a*x for a,x in zip(row,H0vec))%P==0 for row in HA)
    h0_mixed=[]
    for d in data:
        h02=sum(p.get((0,2),0)*x for p,x in zip(d['mons'],H0vec))%P
        if d['kind']=='fonly':
            h11=sum(p.get((1,1),0)*x for p,x in zip(d['mons'],H0vec))%P
            assert h02==0 and h11
            h0_mixed.append(h11)
        if d['kind']=='axis':
            assert sum(p.get((2,1),0)*x for p,x in zip(d['mons'],H0vec))%P==0
    report=dict(certified='exact finite-field minor and nonzero-value witnesses',prime=P,ell=ell,c=c,sqrt_minus_one=ii,
                nodes=[{k:d[k] for k in ['z','w','kind','det']} for d in data],support=support,
                H_rows=len(HA),H=hs,H_corner=corner,H_infinity_linear=H[54],H_infinity_barlinear=H[58],original_J_corner=original_corner,K_corner_prescribable=True,K_rows=len(KA),K=kr,C_rows=len(CA),C=cr,
                axis_inner_coefficients=nonzero_inner,
                minor_determinants=minors,
                H0_coefficients=[[i,j,v] for (i,j),v in sorted(H0.items())],H0_fonly_mixed=h0_mixed,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    target=Path(__file__).with_name('golden_4_11_129_rank_verification.json');target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['prime','ell','c','H_rows','K_rows','C_rows','minor_determinants','H_corner','H_infinity_linear','axis_inner_coefficients','source_sha256']},indent=2),flush=True)

if __name__=='__main__':
    if not __debug__:raise SystemExit('Assertions are required; do not run with -O.')
    main()
