"""Rows of Tables 1-3 of the paper, from data/*.txt and exact arithmetic.
Run from the repository root: python3 code/tables.py"""
import math
from fractions import Fraction
from math import gcd

C, KAPPA = 0.2432599441929454986, 0.7047534517059478841
C2, KAPPA2 = 0.2290907914889382719, 0.0248264240338519476


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def eps(n):
    s = (n & -n).bit_length() - 1
    return {0: 0, 1: -1, 2: 0}.get(s, 1)


def fmt_fact(n):
    return r"\cdot ".join(f"{p}^{{{r}}}" if r > 1 else f"{p}" for p, r in sorted(factor(n).items()))


print("% Table 1")
for n in [2, 6, 10, 12, 24, 26, 35, 45, 98, 105, 120]:
    f = factor(n); om = len(f)
    phi = 1
    for p, r in f.items():
        phi *= p ** (r - 1) * (p - 1)
    sq = len({u * u % n for u in range(1, n) if gcd(u, n) == 1}) if n > 1 else 1
    rho = Fraction(phi, 2 ** om)
    inE = n % 4 == 2 and all(p % 4 == 3 for p in factor(n // 2))
    rs = str(rho.numerator) if rho.denominator == 1 else rf"{rho.numerator}/{rho.denominator}"
    print(f"${n}$ & ${fmt_fact(n)}$ & ${eps(n)}$ & ${phi}$ & ${2**om}$ & ${sq}$ & ${rs}$ & {'yes' if inE else 'no'} \\\\")

print("% Table 2")
for line in open("data/odd_3mod4_counts.txt"):
    y, a = map(int, line.split()); x = 2 * y
    if x in (10**4, 10**5, 10**6, 10**7, 10**8, 10**9, 10**10, 2 * 10**10):
        L = math.log(x)
        one = C * x / math.sqrt(L); two = one * (1 + KAPPA / L)
        print(f"$10^{{{round(math.log10(x))}}}$ & {a} & {one:.0f} & {two:.0f} & {a/one:.5f} & {a/two:.5f} & {(a/two-1)*L*L:.2f} \\\\" if x % 10**round(math.log10(x)) == 0 and str(x).startswith('1') else
              f"$2\\cdot10^{{{round(math.log10(x/2))}}}$ & {a} & {one:.0f} & {two:.0f} & {a/one:.5f} & {a/two:.5f} & {(a/two-1)*L*L:.2f} \\\\")

print("% Table 3")
for line in open("data/rho_sums.txt"):
    x, t = map(int, line.split())
    if x in (10**4, 10**5, 10**6, 10**7, 10**8, 10**9, 10**10):
        S = Fraction(t, 2); L = math.log(x)
        one = C2 * x * x / math.sqrt(L); two = one * (1 + KAPPA2 / L)
        s = f"{float(S):.6e}"
        m, e = s.split("e")
        print(f"$10^{{{round(math.log10(x))}}}$ & ${m}\\cdot10^{{{int(e)}}}$ & {float(S)/one:.6f} & {float(S)/two:.6f} & {(float(S)/two-1)*L*L:.3f} \\\\")

print("% (ratio-1) L^2 ranges for x >= 10^6")
r1 = []
for line in open("data/odd_3mod4_counts.txt"):
    y, a = map(int, line.split()); x = 2 * y
    if x >= 10**6:
        L = math.log(x); r1.append((a / (C * x / math.sqrt(L) * (1 + KAPPA / L)) - 1) * L * L)
print("E:", min(r1), max(r1))
r2 = []
for line in open("data/rho_sums.txt"):
    x, t = map(int, line.split())
    if x >= 10**6:
        L = math.log(x); r2.append((t / 2 / (C2 * x * x / math.sqrt(L) * (1 + KAPPA2 / L)) - 1) * L * L)
print("rho:", min(r2), max(r2))
