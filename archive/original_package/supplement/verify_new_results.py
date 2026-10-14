"""Exact diagnostic checks for the accompanying new polyanalytic results.

Requires Python 3 and SymPy. Run from any directory. The default report is
written next to this script. Finite checks supplement the all-degree proofs.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import math
import sympy as sp

x = sp.Symbol('x')
PRIME = 1000003

def add_polys(B, A, shift, sign):
    C = dict(B)
    for j, c in A.items():
        C[j + shift] = C.get(j + shift, F(0)) + sign * c
    return {j: c for j, c in C.items() if c}

def variations(C):
    signs = [1 if C[j] > 0 else -1 for j in sorted(C)]
    return sum(u != v for u, v in zip(signs, signs[1:]))

def mod_coprime(A, B):
    def modpoly(C):
        vals = {j: int(c.numerator % PRIME) * pow(int(c.denominator % PRIME), -1, PRIME) % PRIME
                for j, c in C.items()}
        p = sp.Poly.from_dict({(j,): v for j, v in vals.items()}, x, modulus=PRIME)
        assert p.degree() == max(C)
        return p
    return sp.gcd(modpoly(A), modpoly(B)).degree() == 0

def dominance(C, delta):
    # Coefficients are based on delta^(exponent^2), with small rational
    # magnitude changes allowed. Test the full inequality exactly.
    for e, c in C.items():
        r = delta ** (-2 * e)
        own = abs(c) * r ** e
        other = sum((abs(v) * r ** j for j, v in C.items() if j != e), F(0))
        assert own > other, (e, 'dominance failed')
    return len(C)

def example(a, b, k):
    assert k + 2*a > 2*b
    delta = F(1, 8 * (k + 2*a + 2*b + 2))
    even_overlap = k % 2 == 0 and k <= 2*b
    if even_overlap:
        h, D = k//2, a+k//2
        Cp = {j: (-1)**j * delta**(j*j) for j in range(D+1)}
        Cm = {j: c * (1 if j < h else -2 if j <= b else -1) for j,c in Cp.items()}
        B = {j: (Cp[j]+Cm[j])/2 for j in range(b+1)}
        A = {j-h: (Cp[j]-Cm[j])/2 for j in range(h,D+1)}
    else:
        exponents = sorted(list(range(0,2*b+1,2)) + list(range(k,k+2*a+1,2)))
        assert len(exponents) == len(set(exponents))
        C = {e: (-1)**j * delta**(e*e) for j,e in enumerate(exponents)}
        B = {j: C[2*j] for j in range(b+1)}
        A = {j: C[k+2*j] for j in range(a+1)}
    if not mod_coprime(A,B):
        old = B[0]
        for j in range(1, a+20):
            B[0] = old + delta**((k+2*a+2*b+4)**2) / j
            if mod_coprime(A,B):
                break
        else:
            raise AssertionError('could not certify coprimality modulo the selected prime')
    if even_overlap:
        plus, minus = add_polys(B,A,k//2,1), add_polys(B,A,k//2,-1)
    else:
        Br, Ar = {2*j:c for j,c in B.items()}, {2*j:c for j,c in A.items()}
        plus, minus = add_polys(Br,Ar,k,1), add_polys(Br,Ar,k,-1)
    comparisons = dominance(plus,delta) + dominance(minus,delta)
    count = variations(plus) + variations(minus)
    target = min(2*(a+b)+1, k+2*a-int(k%2 == 0))
    assert count == target, (a,b,k,count,target)
    assert B[0] and max(A)==a and max(B)==b
    return plus, minus, {'a':a,'b':b,'k':k,'positive_radii':count,
                          'planar_zeros':k*count,'dominance_comparisons':comparisons,
                          'coprimality_prime':PRIME}

def sympoly(C):
    return sp.Poly.from_dict({(j,):sp.Rational(c.numerator,c.denominator) for j,c in C.items()},x)

def main():
    cases = []
    total_comparisons = 0
    for a in range(4):
        for b in range(4):
            for k in range(1,10):
                if k+2*a <= 2*b:
                    continue
                plus,minus,record = example(a,b,k)
                total_comparisons += record['dominance_comparisons']
                cases.append(record)
    sturm = []
    for a,b,k in [(1,2,7),(2,2,3),(2,2,2),(0,2,5)]:
        plus,minus,record = example(a,b,k)
        counts = [int(sympoly(C).count_roots(0,sp.oo)) for C in (plus,minus)]
        assert sum(counts) == record['positive_radii']
        sturm.append({'a':a,'b':b,'k':k,'positive_factor_roots':counts})
    grid = 0
    for m in range(1,101):
        for n in range(m+1,10*m+11):
            grid += 1
            U = min(m*m+n*n,(4*m+1)*n-5*m)
            V = min(n*n,(4*m+1)*n-2*m*(m+2))
            if n >= 2*m+1:
                s = max(0,3*m+1-n)
                L = (4*m+1)*(n-m)+max(0,3*m-n)*max(0,3*m+1-n)
                assert L == (n-m+s)*(4*m-2*s+1)+s*s
                assert U-4*m*(m-1) <= L <= U
            else:
                assert U == m*m+n*n
            if n <= 2*m+1:
                E_lower = n*n
            else:
                a = min(m,(n-2*m)//4)
                k = n-2*a-1
                assert k > 2*m
                constructed = n+2*(m+a)*(n-2*a-1)
                admissible = [n+2*(m+j)*(n-2*j-1) for j in range(m+1) if n-2*j-1>2*m]
                assert constructed == max(admissible)
                # Also compare against every overlapping-support version of
                # this particular completion template.
                all_templates = [(n-2*j-1)*min(2*(j+m)+1,n-1-int(n%2==1))+2*j+1 for j in range(m+1)]
                assert constructed == max(all_templates)
                E_lower = max((2*m+1)*n,constructed)
            assert V-6*m*m <= E_lower <= V
            if n==4*m+2:
                assert E_lower == 9*m*m+7*m+2
                assert E_lower-((4*m+1)*n-4*m*(2*m+1)) == m*(m-1)
    report = {'status':'passed', 'arithmetic':'exact rational and integer arithmetic',
              'radial_cases':len(cases), 'dominance_comparisons':total_comparisons,
              'coprimality_certificates':len(cases), 'sturm_cases':sturm,
              'integer_degree_pairs':grid, 'integer_range':'1 <= m <= 100; m < n <= 10m+10',
              'cases':cases,
              'scope':'Checks new formulas only. Infinite-degree conclusions use the accompanying proofs. Original certificate files were not supplied or rerun.'}
    path = Path(__file__).with_name('new_results_verification.json')
    path.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='cases'},indent=2))

if __name__ == '__main__':
    main()
