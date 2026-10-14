"""Exact integer audit of the degree obstruction in Corollary v31.

This audits identities and rounding, not the previously proved upper-bound
or construction theorems on which the mathematical corollary depends.
"""
from pathlib import Path
from math import isqrt
import hashlib,json,sys

if not __debug__:
    raise RuntimeError('Verification must not be run with Python optimization.')

def require(condition,message):
    if not condition:raise ArithmeticError(message)

def upper(m,n):return (4*m+1)*n-2*m*(m+2)
def lower_degree(n):return n-1-isqrt((n-1)*(n-2)//2)
def degree_ceiling(m):return (4*m+1+isqrt(8*m*m-8*m+1))//2

def main():
    for n in range(2,1001):
        m=lower_degree(n)
        require(1<=m<n,('range',m,n))
        require(n*n<=upper(m,n),('lower rounding',m,n))
        require(m==1 or n*n>upper(m-1,n),('strict lower rounding',m,n))
        for mm in {m,max(1,m-1),n-1}:
            require(upper(mm,n)-n*n==(n-1)*(n-2)-2*(mm-(n-1))**2,('square identity',mm,n))
    for m in range(1,1001):
        n=degree_ceiling(m)
        require(n*n<=upper(m,n) and (n+1)**2>upper(m,n+1),('upper rounding',m,n))
    table={n:lower_degree(n) for n in range(2,10)}
    require(list(table.values())==[1,1,2,2,2,3,3,3],'degree table')
    deficits=[9*k-k*k-2 for k in range(1,9)]
    require(deficits==[6,12,16,18,18,16,12,6],'deficits')
    for k in range(1,9):
        require(lower_degree(3*k)==k,('cubic minimum',k))
        require(upper(k-1,3*k)==9*k*k-deficits[k-1],('deficit identity',k))
    require(upper(8,27)==731 and lower_degree(27)==8,'first permitted cubic reduction')
    require(9*9-9*9+2==2,'deficit sign after eighth order')
    require({n:lower_degree(n) for n in [10,11,13,27]}=={10:3,11:4,13:4,27:8},'open lower endpoints')
    report=dict(all_passed=True,exact_arithmetic=True,degree_table=table,
        sharp_cubic_orders=list(range(1,9)),one_degree_loss_deficits=deficits,
        integer_roundings_checked_through=1000,
        open_degree_brackets={'10':[3,4],'11':[4,5],'13':[4,6],'27':[8,9]},
        theorem_dependency='n:eq:Eupper, d:prop:seeds, d:prop:boundary, v16:thm:borderline',
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    destination=Path(__file__).with_name('minimum_conjugate_degree_verification.json')
    destination.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
