"""Exact standard-library verification of the printed P40 and P49.

Default: verify both adjacent sharp_2_6_40_certificate.json and
sharp_2_7_49_certificate.json. Optional positional paths override them.
Every proof inequality uses fractions.Fraction; floats are display only.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,time,hashlib

if not __debug__:
    raise RuntimeError('Run without -O: all exact assertions must be enabled.')

DATA={
40:(10**10,
[[9989148,-1437879431,42110613567],[-1488721221,85880077969,44531370007],
 [41495821383,181048684005,-8945073549],[127109446221,24463634439,-54584415742],
 [70927784116,-90270783729,-13319637790],[-30642879880,-58978230217,12977971008],
 [171904696,8861758778,10000000000]],
[[2386743,-536302879,19849756729],[-170102520,20587381168,26848909014],
 [1758707049,27892403146,9999884021],[-18670288098,19443937659,-6680374702],
 [-10494108981,26168926864,-13516927630],[7007756920,1936949230,-10640650833],
 [4236054207,1869520737,0]]),
49:(1000,
[[-42,-176,-139],[-227,933,-432],[-179,-1457,737],[794,3712,353],
 [-2958,3048,-1571],[7455,1348,-1396],[-448,-5235,-166],[-68,158,1000]],
[[-33,-479,-2287],[-262,-1987,653],[-1501,-8252,965],[-4456,4656,3427],
 [-7424,1665,-1101],[1057,8157,-195],[-1669,-264,-2045],[363,697,0]])}

def add(a,b):return(a[0]+b[0],a[1]+b[1])
def mul(a,b):return(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(a,b):return(a*b[0],a*b[1])
def conj(a):return(a[0],-a[1])

def evaluate(coeffs,z,w):
    zp=[(F(1),F(0))];wp=zp.copy()
    for _ in range(max(i for i,j,c in coeffs)):zp.append(mul(zp[-1],z))
    for _ in range(max(j for i,j,c in coeffs)):wp.append(mul(wp[-1],w))
    f=az=aw=(F(0),F(0))
    for i,j,c in coeffs:
        f=add(f,mul(c,mul(zp[i],wp[j])))
        if i:az=add(az,scale(i,mul(c,mul(zp[i-1],wp[j]))))
        if j:aw=add(aw,scale(j,mul(c,mul(zp[i],wp[j-1]))))
    return f,az,aw

def inverse(matrix):
    n=len(matrix)
    a=[list(row)+[F(int(i==j)) for j in range(n)] for i,row in enumerate(matrix)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if a[i][j]),None)
        assert pivot is not None,'Singular center Jacobian'
        a[j],a[pivot]=a[pivot],a[j]
        p=a[j][j];a[j]=[v/p for v in a[j]]
        for i in range(n):
            if i!=j:
                p=a[i][j];a[i]=[u-p*v for u,v in zip(a[i],a[j])]
    return [row[n:] for row in a]

def trim(p):
    p=[int(a)%101 for a in p]
    while len(p)>1 and p[-1]==0:p.pop()
    return p

def pgcd(a,b):
    a,b=trim(a),trim(b)
    while b!=[0]:
        r=a.copy()
        while len(r)>=len(b) and r!=[0]:
            d=len(r)-len(b);v=r[-1]*pow(b[-1],-1,101)%101
            for j,c in enumerate(b):r[j+d]=(r[j+d]-v*c)%101
            r=trim(r)
        a,b=b,r
    return [v*pow(a[-1],-1,101)%101 for v in a]

def verify(path):
    raw=json.loads(path.read_text(encoding='utf-8'))
    N=int(raw['polynomial'][1:]);assert N in DATA
    cd,real,imag=DATA[N];n=len(real)-1
    assert raw['real_numerators']==real and raw['imaginary_numerators']==imag,'Printed coefficient mismatch'
    assert int(raw['coefficient_denominator'])==cd
    assert (raw['analytic_degree'],raw['conjugate_degree'])==(n,2)
    coeffs=[(i,j,(F(real[i][j],cd),F(imag[i][j],cd))) for i in range(n+1) for j in range(3)]
    assert coeffs[-1][2]==(F(1),F(0))
    star=[(j,i,conj(c)) for i,j,c in coeffs]
    p0=[real[i][0]+10*imag[i][0] for i in range(n+1)]
    p1=[real[i][1]+10*imag[i][1] for i in range(n+1)]
    u=[real[0][j]+10*imag[0][j] for j in range(3)]
    v=[real[i][0]-10*imag[i][0] for i in range(n+1)]
    assert all(a[-1]%101 for a in [p0,p1,u,v]),'Degree drop modulo 101'
    assert pgcd(p0,p1)==[1] and pgcd(u,v)==[1],'Coprimality certificate failed'
    den=int(raw['center_denominator']);assert den==10**40
    r=F(int(raw['radius'][0]),int(raw['radius'][1]));assert r==F(1,10**20)
    rc=[tuple(F(int(v),den) for v in row) for row in raw['real_centers']]
    nc=[tuple(F(int(v),den) for v in row) for row in raw['nonreal_centers']]
    assert len(rc)==N and all(len(v)==2 for v in rc)
    assert len(nc)==n*n+4-N and all(len(v)==4 for v in nc)
    embedded=[(x,y,x,-y) for x,y in rc]+nc
    for i,a in enumerate(embedded):
        for b in embedded[i+1:]:
            assert max(abs(x-y) for x,y in zip(a,b))>2*r,'Overlapping intersection boxes'
    for a,b,c,d in nc:
        assert max(abs(a-c),abs(b+d))>2*r,'Nonreal box meets real slice'
    qmax=brmax=F(0);signs=[];start=time.time()
    for index,center in enumerate(rc+nc):
        isreal=len(center)==2
        z=center[:2];w=conj(z) if isreal else center[2:]
        f,A,B=evaluate(coeffs,z,w)
        if isreal:
            matrix=[[A[0]+B[0],B[1]-A[1]],[A[1]+B[1],A[0]-B[0]]]
            value=list(f)
            det=matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
            assert det
            signs.append(1 if det>0 else -1)
        else:
            g,C,D=evaluate(star,z,w)
            matrix=[[A[0],-A[1],B[0],-B[1]],[A[1],A[0],B[1],B[0]],
                    [C[0],-C[1],D[0],-D[1]],[C[1],C[0],D[1],D[0]]]
            value=list(f)+list(g)
        Y=inverse(matrix)
        norm=max(sum(abs(v) for v in row) for row in Y)
        beta=max(abs(sum(a*b for a,b in zip(row,value))) for row in Y)
        R=max(sum(abs(v) for v in z),sum(abs(v) for v in w))+2*r
        H=sum((i+j)*(i+j-1)*(abs(c[0])+abs(c[1]))*R**(i+j-2)
              for i,j,c in coeffs if i+j>=2)
        q=(4 if isreal else 16)*r*norm*H
        assert q+beta/r<F(1,10**7),('Contraction/self-mapping failed',index)
        qmax=max(qmax,q);brmax=max(brmax,beta/r)
    assert (signs.count(1),signs.count(-1))==((22,18) if N==40 else (27,22))
    report={'polynomial':f'P{N}','all_checks_passed':True,
            'exact_real_slice_zeros':N,'positive':signs.count(1),'negative':signs.count(-1),
            'nonreal_complex_intersections':len(nc),'bezout_exhaustion':n*n+4,
            'coprimality_modulo_101':True,'printed_coefficients_match':True,
            'strict_exact_bound':'q + beta/r < 10^-7 in every box',
            'display_only_q_max':float(qmax),'display_only_beta_over_r_max':float(brmax),
            'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'seconds':round(time.time()-start,3)}
    print(json.dumps(report,indent=2),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('certificates',nargs='*',type=Path)
    p.add_argument('--report',type=Path)
    a=p.parse_args()
    paths=a.certificates or [Path(__file__).with_name(f'sharp_2_{n}_{N}_certificate.json') for n,N in [(6,40),(7,49)]]
    report={'arithmetic':'fractions.Fraction and finite fields; exact throughout',
            'certificates':[verify(path) for path in paths],
            'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if a.report:a.report.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
