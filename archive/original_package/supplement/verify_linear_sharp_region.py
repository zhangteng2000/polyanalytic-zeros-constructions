"""Exact finite diagnostics for the simplified splitting and linear sharp cone.

Run with Python 3 and SymPy, without -O. These checks supplement the proofs
in polyanalytic_sharp_v3.tex. The specialization epsilon=1 checks polynomial
support only; it is NOT an extremizing parameter or a numerical root count.
"""
from pathlib import Path
import json
import time
import sympy as S

z, w, t = S.symbols('z w t')


def global_check(k, r):
    start = time.time()
    N = 2*r + 1
    p = lambda e: S.Poly(e, z, w, t, domain=S.QQ_I)
    lam = 1 + 2*S.I
    Q = p(1)
    for h in range(1, r+1):
        Q *= p(z*w - h*h)
    f = p(z**N) + Q
    fs = p(w**N) + Q
    B = p(-lam*z**r)
    Bs = p(-S.conjugate(lam)*w**r)
    A = p(lam*w**(r+1)) - B*p(z*w)
    G = A*f + B*fs
    assert G.degree(w) <= 2*r
    assert G.degree(z) <= 3*r+2
    assert max(i+j for (i,j,h), c in G.terms()) <= 3*r+3
    betaV = p(t*S.Rational(1,k))*G
    gamma = p(t**3*S.Rational(1,k**3))*B*B*Bs
    L = [p(3), f, f*f + 2*betaV]
    for h in range(3, k+1):
        L.append(f*L[-1] + betaV*L[-2] - gamma*L[-3])
    P = L[k] + (-p(t*S.Rational(1,k))*Bs)**k
    n, m = k*N, k*r
    assert P.degree(w) <= m and P.degree(z) == n
    assert max(i+j for (i,j,h), c in P.terms()) == n
    assert P.coeff_monomial(z**n) == 1
    assert P.coeff_monomial(z**m*w**m) == 1
    delta = P - f**k
    assert all(h >= 1 for (i,j,h), c in delta.terms())
    if r >= 2:
        assert all(i+j < n for (i,j,h), c in delta.terms())
    else:
        # Near-pure is all that is claimed at r=1; mixed highest terms occur.
        assert any(i+j == n and h > 0 and j > 0
                   for (i,j,h), c in delta.terms())
    return {'k': k, 'r': r, 'degree_z': int(P.degree(z)),
            'degree_w': int(P.degree(w)), 'total_degree': n,
            'pure_head': r >= 2, 'terms': len(P.terms()),
            'seconds': round(time.time()-start, 3)}


def local_scaling_check(k):
    u, v, b, c, eta, tau = S.symbols('u v b c eta tau')
    p = lambda e: S.Poly(e, u, v, b, c, eta, tau, domain=S.QQ)
    U, B, C, T = p(u), p(b), p(c), p(tau)
    V = p(v+eta*u)
    betaV = B*V*S.Rational(1,k)
    gamma = C*B*B*S.Rational(1,k**3)
    model = [p(3), U, U*U+2*betaV]
    global_local = [p(3), T*U, (T*U)**2+2*T*T*betaV]
    for h in range(3, k+1):
        model.append(U*model[-1]+betaV*model[-2]-gamma*model[-3])
        global_local.append(T*U*global_local[-1]+T*T*betaV*global_local[-2]
                            -T**3*gamma*global_local[-3])
    H = model[k] + (-C*S.Rational(1,k))**k
    P = global_local[k] + (-T*C*S.Rational(1,k))**k
    assert (P-T**k*H).is_zero
    assert all(mon[0]+mon[1] < k for mon, coeff in (H-U**k).terms())
    assert all(mon[0]+mon[1]+mon[2]+mon[3] == k
               for mon, coeff in H.terms())


def shear_identity_check():
    # The identity holds modulo the seed ideal for symbolic radial coefficients.
    reports = []
    for r in range(1, 9):
        N = 2*r+1
        lam = 1+2*S.I
        Q = S.prod(z*w-S.Rational(j*j, 10000) for j in range(1,r+1))
        f = z**N+Q
        fs = w**N+Q
        B = -lam*z**r
        A = lam*w**(r+1)-B*z*w
        # (zw)^r (A/B+zw) = -w^N = Q-f*.
        lhs = S.cancel((z*w)**r*(A/B+z*w))
        assert S.expand(lhs-Q+fs) == 0
        reports.append(r)
    return reports


def inversion_support_check():
    tested = 0
    for n in range(2, 31):
        for m in range(1, n):
            for j in range(m+1):
                for i in range(n-j+1):
                    # z^i w^j -> u^(n-i) v^(m-j), normalized by s^m.
                    assert n-i >= 0 and m-j >= 0
                    scale_power = n-i-j
                    assert scale_power >= 0
                    assert (scale_power == 0) == (i+j == n)
                    if scale_power == 0:
                        assert n-i == j
                    tested += 1
    return tested


def cone_check():
    tested = 0
    for m in range(1, 1001):
        for n in range(m+1, 3*m+1):
            if n == 3*m-1:
                continue
            s = n % 3
            k = (n-s)//3
            a = m-k-s
            assert k >= 1 and a >= 0 and 0 <= s <= 2
            assert k+a < 3*k
            assert (k+a+s, 3*k+s) == (m,n)
            assert 3*m-n == 3*a+2*s
            base = (k+a)**2+(3*k)**2
            assert base+2*s*(4*k+a+s) == m*m+n*n
            tested += 1
    return tested


def run():
    start = time.time()
    local_orders = list(range(2, 21))
    for k in local_orders:
        local_scaling_check(k)
    print('Local scaling and sheared pure-head identities passed.', flush=True)
    cases = [(k,r) for r in range(1,5) for k in range(2,9)]
    cases += [(12,1),(12,2)]
    degrees = []
    for k,r in cases:
        result = global_check(k,r)
        degrees.append(result)
        print(json.dumps(result), flush=True)
    shear = shear_identity_check()
    inversion = inversion_support_check()
    cone = cone_check()
    report = {
        'status': 'passed',
        'scope': 'Finite exact diagnostics, not a universal proof or certified parameter choice.',
        'arithmetic': 'Exact integers, rationals, and Gaussian rationals; no numerical roots.',
        'local_model_orders': local_orders,
        'global_degree_cases': degrees,
        'shear_identity_r': shear,
        'inversion_monomials': inversion,
        'cone_pairs': cone,
        'cone_max_m': 1000,
        'excluded_cone_line': 'n=3m-1; no impossibility claim is made.',
        'seconds': round(time.time()-start, 3)
    }
    target = Path(__file__).with_name('linear_sharp_region_verification.json')
    target.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'passed', 'degree_cases': len(degrees),
                      'cone_pairs': cone, 'inversion_monomials': inversion,
                      'seconds': report['seconds']}), flush=True)


if __name__ == '__main__':
    run()
