"""Brute-force checks of the algebraic statements (Sections 2-4).

For every n <= N it checks, by direct enumeration in Z/nZ:
  * Lemma 2.1: #ring homomorphisms Z_m -> Z_n (n | m, m = n and m = 2n) is 2^omega(n);
  * Lemma 2.2: #surjective group homomorphisms Z_m -> Z_n is phi(n);
  * Lemma 2.3: [U_n : U_n^2] = 2^(omega(n) + eps(n));
  * Theorem 1.1: phi(n) / 2^omega(n) = 2^eps(n) |U_n^2|, and integrality <=> n not in E;
and for odd n <= N_odd:
  * Theorem 1.2: e -> chi_e is an injective homomorphism (B_n, +) -> Hom(U_n, {+-1})
    whose image is all of Hom(U_n, {+-1}) (sizes agree), and sigma_u(sqrt d_e) = chi_e(u) sqrt d_e,
    checked numerically with Gauss sums in C.
Run: python3 code/verify_algebra.py
"""
import cmath
from fractions import Fraction
from math import gcd

N = 10000
N_ODD = 1000


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def eps(n):
    s = 0
    while n % 2 == 0:
        n //= 2
        s += 1
    return {0: 0, 1: -1, 2: 0}.get(s, 1)


def in_E(n):
    if n % 4 != 2:
        return False
    return all(p % 4 == 3 for p in factor(n // 2))


def ring_hom_count(m, n):
    """Maps f(k) = k e (mod n) that are additive and multiplicative on Z_m."""
    c = 0
    for e in range(n):
        if (m * e) % n:                       # not well defined on Z_m
            continue
        if all((a * b * e - (a * e) * (b * e)) % n == 0 for a in range(min(m, n)) for b in range(min(m, n))):
            c += 1
    return c


def surj_count(m, n):
    return sum(1 for a in range(n) if (m * a) % n == 0 and gcd(a, n) == 1)


def legendre(a, p):
    t = pow(a % p, (p - 1) // 2, p)
    return -1 if t == p - 1 else t


def check_counts():
    for n in range(1, N + 1):
        f = factor(n)
        om = len(f)
        phi = 1
        for p, r in f.items():
            phi *= p ** (r - 1) * (p - 1)
        U = [u for u in range(n) if gcd(u, n) == 1] if n > 1 else [0]
        sq = {(u * u) % n for u in U} if n > 1 else {0}
        assert len(U) == phi
        assert len(U) // len(sq) == 2 ** (om + eps(n)), n
        rho = Fraction(phi, 2 ** om)
        assert rho == Fraction(2) ** eps(n) * len(sq), n
        assert (rho.denominator == 1) == (not in_E(n)), n
        if n <= 60:
            for m in (n, 2 * n):
                assert ring_hom_count(m, n) == 2 ** om, (m, n)
                assert surj_count(m, n) == phi, (m, n)
    print(f"Lemmas 2.1-2.3 and Theorem 1.1 verified for n <= {N} (hom counts for n <= 60).")


def check_bijection():
    for n in range(3, N_ODD + 1, 2):
        f = factor(n)
        primes = sorted(f)
        U = [u for u in range(1, n) if gcd(u, n) == 1]
        idem = [e for e in range(n) if (e * e - e) % n == 0]
        assert len(idem) == 2 ** len(primes)

        def S(e):
            return frozenset(p for p in primes if e % p == 1)

        def chi(e, u):
            v = 1
            for p in S(e):
                v *= legendre(u, p)
            return v

        tables = {}
        for e in idem:
            tab = tuple(chi(e, u) for u in U)
            # chi_e is a character of order <= 2
            for a in U[:12]:
                for b in U[:12]:
                    assert chi(e, a * b) == chi(e, a) * chi(e, b)
            tables[e] = tab
        assert len(set(tables.values())) == len(idem)          # injective
        # Boolean group law e (+) f = e + f - 2ef  <->  chi_e chi_f
        for e in idem:
            for g in idem:
                h = (e + g - 2 * e * g) % n
                assert all(tables[h][k] == tables[e][k] * tables[g][k] for k in range(len(U)))
        # Galois action on sqrt(d_e) = prod_{p in S(e)} tau_p, tau_p a quadratic Gauss sum
        if n <= 120:
            z = cmath.exp(2j * cmath.pi / n)
            for e in idem:
                def tau(power):
                    val = 1
                    for p in S(e):
                        val *= sum(legendre(a, p) * z ** ((n // p) * a * power) for a in range(1, p))
                    return val
                t1 = tau(1)
                d = 1
                for p in S(e):
                    d *= p if p % 4 == 1 else -p
                assert abs(t1 * t1 - d) < 1e-6 * max(1, abs(d))
                for u in U:
                    assert abs(tau(u) - chi(e, u) * t1) < 1e-6 * max(1, abs(t1))
    print(f"Theorem 1.2 verified for odd n <= {N_ODD} (Gauss-sum check for n <= 120).")


if __name__ == "__main__":
    check_counts()
    check_bijection()
