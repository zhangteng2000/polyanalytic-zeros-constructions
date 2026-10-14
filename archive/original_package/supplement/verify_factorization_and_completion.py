"""Finite exact-arithmetic diagnostics for the v8 analytic proofs.

These checks do not prove the all-degree topological constructions or
identify any of the unresolved irreducible maxima. They independently
enumerate factor-degree bounds and optimize the stated seed set.
"""
from functools import lru_cache
from pathlib import Path
import json
import sympy as sp


def upper(a, b):
    a, b = sorted((a, b))
    if a == 0:
        return b
    if a == b:
        return 2 * a * a
    return min(a * a + b * b, (4 * a + 1) * b - 5 * a)


def seed_admissible(a, b, m, n):
    return (1 <= a <= m and a < b <= 3 * a and a + b < n
            and not (a >= 4 and b == 3 * a - 1))


def seed_count(a, b, n):
    return n + a * a + b * b + a - b


def candidates(m, n):
    a0 = min(m, (n - 1) // 4)
    answer = []
    for a in range(a0, min(m, a0 + 2) + 1):
        b = min(3 * a, n - a - 1)
        if a >= 4 and b == 3 * a - 1:
            b -= 1
        if a >= 1 and b > a:
            assert seed_admissible(a, b, m, n)
            answer.append((seed_count(a, b, n), a, b))
    return answer


def radial(m, n):
    a = min(m, (n - 2 * m) // 4)
    return max((2 * m + 1) * n,
               n + 2 * (m + a) * (n - 2 * a - 1))


def main():
    report = {"status": "running", "scope": __doc__.strip()}
    # This evaluates both individual upper bounds, not a copy of the
    # concavity proof used in the paper.
    split_checks = 0
    for m in range(2, 32):
        for n in range(m + 1, 101):
            target = upper(m - 1, n) + 1
            for a in range(1, m):
                for b in range(n + 1):
                    assert upper(a, b) + upper(m - a, n - b) <= target
                    split_checks += 1
    report["split_conjugate_degree_inequalities"] = split_checks

    a, b, m, n, s = sp.symbols("a b m n s", integer=True)
    symbols = [
        ((m - 1) ** 2 + 1 - (a + (m - a) ** 2)
         - (a - 1) * (2 * m - a - 2)),
        ((4 * m - 3) * n - 5 * m + 6
         - (4 * a * b + 4 * (m - a) * (n - b) + m + n - 6 * s)
         - (4 * ((a - 1) * n + (m - 2 * a) * b) + 6 * (s - m + 1))),
        ((4 * m + 1) * (n - m - 1) + 1
         - ((4 * m - 3) * n - 5 * m + 6)
         - (4 * n - 4 * m * m - 6)),
        ((a + 1) ** 2 + (n - a - 2) ** 2 + 2 * (a + 1) + 1
         - (a * a + (n - a - 1) ** 2 + 2 * a + 1)
         - (4 * a - 2 * n + 6)),
    ]
    assert all(sp.expand(value) == 0 for value in symbols)
    report["symbolic_identities"] = len(symbols)

    envelope_checks = 0
    seed_checks = 0
    for m in range(1, 33):
        for n in range(max(5, m + 1), 151):
            brute = []
            for a in range(1, m + 1):
                for b in range(a + 1, 3 * a + 1):
                    if seed_admissible(a, b, m, n):
                        brute.append(seed_count(a, b, n))
                        seed_checks += 1
            calculated = max(value for value, _, _ in candidates(m, n))
            assert calculated == max(brute), (m, n, calculated, max(brute))
            assert calculated <= min(n * n, (4 * m + 1) * n - 2 * m * (m + 2))
            envelope_checks += 1
    report["completion_envelopes_enumerated"] = envelope_checks
    report["individual_seed_transfers_enumerated"] = seed_checks

    subdiagonals = 0
    gains = 0
    for m in range(4, 1001):
        formulas = [(4 * m - 2, 10 * m * m - 16 * m + 10),
                    (4 * m - 1, 10 * m * m - 10 * m + 5),
                    (4 * m, 10 * m * m - 10 * m + 6)]
        for n, value in formulas:
            assert max(t[0] for t in candidates(m, n)) == value
            subdiagonals += 1
        if m >= 11:
            expected_gains = [(m - 1) * (m - 10),
                              m * m - 8 * m + (6 if m % 2 == 0 else 5),
                              m * m - 11 * m + 6]
            for (n, value), expected in zip(formulas, expected_gains):
                assert value - radial(m, n) == expected > 0
                gains += 1
        n = 4 * m + 1
        assert max(t[0] for t in candidates(m, n)) == n + 10 * m * m - 2 * m
    report["subdiagonal_formulas"] = subdiagonals
    report["strict_radial_improvement_checks"] = gains

    @lru_cache(None)
    def max_sum_squares(total, count):
        if count == 1:
            return total * total
        return max(d * d + max_sum_squares(total - d, count - 1)
                   for d in range(1, total - count + 2))

    factor_checks = 0
    for n in range(2, 101):
        for ell in range(1, n + 1):
            assert max_sum_squares(n, ell) == (n - ell + 1) ** 2 + ell - 1
            factor_checks += 1
    report["factor_degree_partitions_optimized"] = factor_checks

    report["examples"] = [
        {"m": 100, "n": n, "radial": radial(100, n),
         "completion": max(t[0] for t in candidates(100, n))}
        for n in (398, 399, 400)
    ]
    assert [r["completion"] for r in report["examples"]] == [98410, 99005, 99006]
    assert [r["radial"] for r in report["examples"]] == [89500, 89799, 90100]
    report["status"] = "passed"
    output = Path(__file__).with_name("factorization_and_completion_verification.json")
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
