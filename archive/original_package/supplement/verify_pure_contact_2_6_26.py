"""Exact rational verification of the 26-zero pure contact-six seed.
No numerical root finder is used. The coefficients in the certificate
define the polynomial; every required support and contact condition is checked.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,time
if not __debug__:raise RuntimeError('Do not run with Python -O.')
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def mul(a,b):return(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(a,b):return(a*b[0],a*b[1])
def conj(a):return(a[0],-a[1])
def inverse(matrix):
    n=len(matrix);a=[list(row)+[F(int(i==j)) for j in range(n)] for i,row in enumerate(matrix)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if a[i][j]),None);assert pivot is not None
        a[j],a[pivot]=a[pivot],a[j];v=a[j][j];a[j]=[x/v for x in a[j]]
        for i in range(n):
            if i!=j:
                v=a[i][j];a[i]=[x-v*y for x,y in zip(a[i],a[j])]
    return [row[n:] for row in a]
def evaluate(coeffs,z,w):
    zp=[(F(1),F(0))];wp=zp.copy()
    for _ in range(6):zp.append(mul(zp[-1],z))
    for _ in range(6):wp.append(mul(wp[-1],w))
    f=A=B=(F(0),F(0))
    for i,j,c in coeffs:
        f=add(f,mul(c,mul(zp[i],wp[j])))
        if i:A=add(A,scale(i,mul(c,mul(zp[i-1],wp[j]))))
        if j:B=add(B,scale(j,mul(c,mul(zp[i],wp[j-1]))))
    return f,A,B
EXPECTED_COEFFICIENTS={(0, 1): ('1095184937584699/1125899906842624', '-4000422529849379/18014398509481984'), (0, 2): ('3530749276274985/18014398509481984', '-373254250500697/9007199254740992'), (1, 1): ('2635113859063031/4503599627370496', '-5868212651896887/9007199254740992'), (1, 2): ('-6872376342183993/2251799813685248', '2500801915392709/2251799813685248'), (2, 1): ('-6173249439061019/4503599627370496', '-1884493385596669/9007199254740992'), (2, 2): ('-1402124144343637/1125899906842624', '8597653020315327/4503599627370496'), (3, 1): ('-151946931462199/562949953421312', '-3139336158634565/4503599627370496'), (3, 2): ('4669734454329891/1125899906842624', '3680316199713847/4503599627370496'), (4, 1): ('-3283453096864763/2251799813685248', '-8448110744614189/72057594037927936'), (6, 0): ('1', '0')}

def verify(path,min_count):
    raw=json.loads(path.read_text(encoding='utf-8'));assert raw['kind']=='pure_contact_2_6_lower_bound'
    coeffs=[(int(i),int(j),(F(a),F(b))) for i,j,a,b in raw['terms']]
    assert len(set((i,j) for i,j,c in coeffs))==len(coeffs)
    assert {(i,j):(str(c[0]),str(c[1])) for i,j,c in coeffs}==EXPECTED_COEFFICIENTS,'Printed coefficient mismatch'
    assert 'nonreal_centers' in raw and len(raw['real_centers'])==26
    assert all(c!=(0,0) for i,j,c in coeffs)
    assert [(i,j,c) for i,j,c in coeffs if j==0]==[(6,0,(F(1),F(0)))], 'Wrong pure analytic jet'
    assert all(0<=i<=6 and 0<=j<=2 and (i+j<=5 or (i,j)==(6,0)) for i,j,c in coeffs)
    assert any(j==2 for i,j,c in coeffs),'Conjugate degree must be exactly two'
    assert any((i,j)==(0,1) for i,j,c in coeffs),'Origin is not regular'
    den=int(raw['center_denominator']);assert den==10**180
    radius=F(int(raw['radius'][0]),int(raw['radius'][1]));assert radius==F(1,10**80)
    centers=[tuple(F(int(a),den) for a in row) for row in raw['real_centers']]
    assert len(centers)==raw['proposed_count'] and min_count<=len(centers)<=34
    assert all(len(v)==2 for v in centers)
    assert (F(0),F(0)) in centers,'The contact zero must be included'
    nc=[tuple(F(int(a),den) for a in row) for row in raw.get('nonreal_centers',[])]
    if 'nonreal_centers' in raw:
        assert len(nc)+len(centers)==36 and all(len(v)==4 for v in nc)
        for a,b,c,d in nc:assert max(abs(a-c),abs(b+d))>2*radius,'Nonreal box meets real slice'
    allboxes=[(x,y,x,-y) for x,y in centers]+nc
    for k,a in enumerate(allboxes):
        for b in allboxes[k+1:]:assert max(abs(a[j]-b[j]) for j in range(4))>2*radius,'Overlapping boxes'
    signs=[];qmax=betamax=F(0)
    for h,z in enumerate(centers):
        f,A,B=evaluate(coeffs,z,conj(z))
        M=[[A[0]+B[0],B[1]-A[1]],[A[1]+B[1],A[0]-B[0]]]
        det=M[0][0]*M[1][1]-M[0][1]*M[1][0];assert det,'Singular center Jacobian'
        Y=[[M[1][1]/det,-M[0][1]/det],[-M[1][0]/det,M[0][0]/det]]
        norm=max(sum(abs(a) for a in row) for row in Y)
        beta=max(abs(sum(a*b for a,b in zip(row,f))) for row in Y)
        R=abs(z[0])+abs(z[1])+2*radius
        H=sum((i+j)*(i+j-1)*(abs(c[0])+abs(c[1]))*R**(i+j-2) for i,j,c in coeffs if i+j>=2)
        q=4*radius*norm*H
        assert q+beta/radius<F(1,1024),('Contraction failed',h)
        signs.append(1 if det>0 else -1);qmax=max(qmax,q);betamax=max(betamax,beta/radius)
    star=[(j,i,conj(c)) for i,j,c in coeffs]
    for h,v in enumerate(nc):
        z=v[:2];w=v[2:];f,A,B=evaluate(coeffs,z,w);g,C,D=evaluate(star,z,w)
        M=[[A[0],-A[1],B[0],-B[1]],[A[1],A[0],B[1],B[0]],
           [C[0],-C[1],D[0],-D[1]],[C[1],C[0],D[1],D[0]]]
        Y=inverse(M);norm=max(sum(abs(x) for x in row) for row in Y)
        beta=max(abs(sum(a*b for a,b in zip(row,list(f)+list(g)))) for row in Y)
        R=max(sum(abs(x) for x in z),sum(abs(x) for x in w))+2*radius
        H=sum((i+j)*(i+j-1)*(abs(c[0])+abs(c[1]))*R**(i+j-2) for i,j,c in coeffs if i+j>=2)
        q=16*radius*norm*H
        assert q+beta/radius<F(1,1024),('Complex contraction failed',h)
        qmax=max(qmax,q);betamax=max(betamax,beta/radius)
    report={'all_exact_checks_passed':True,'certified_regular_zero_lower_bound':len(centers),
            'positive_boxes':signs.count(1),'negative_boxes':signs.count(-1),
            'contact_order':6,'origin_regular':True,'pure_head':'z^6','exact_conjugate_degree':2,
            'strict_bound':'q+beta/r < 1/1024 in every rational box',
            'display_only_q_max':float(qmax),'display_only_beta_over_r_max':float(betamax),
            'consequent_pure_E_2_7_lower_bound':len(centers)+13,
            'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    if 'nonreal_centers' in raw:
        report.update(exact_real_slice_zeros=len(centers),nonreal_complex_intersections=len(nc),pure_bezout_exhaustion=36)
    return report
def main():
    p=argparse.ArgumentParser();p.add_argument('certificate',type=Path,nargs='?',default=Path(__file__).with_name('pure_contact_2_6_26_certificate.json'));p.add_argument('--min-count',type=int,default=26);p.add_argument('--report',type=Path,default=Path(__file__).with_name('pure_contact_2_6_26_verification.json'))
    a=p.parse_args();r=verify(a.certificate,a.min_count);r['verifier_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if a.report:a.report.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(r,indent=2))
if __name__=='__main__':main()
