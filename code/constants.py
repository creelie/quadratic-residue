"""High-precision values of the constants in Theorems 1.3 and 1.4 (prime zeta
function expansions; no truncation of Euler products).  Run: python3 code/constants.py"""
import mpmath as mp

mp.mp.dps = 40
chi4 = [0, 1, 0, -1]


def R(s):
    return (1 - mp.power(2, -s)) * mp.zeta(s) / mp.dirichlet(s, chi4)


def logP3(s):
    """log prod_{p=3 mod 4}(1-p^-s)^-1 = sum_k 2^-k log R(2^(k-1) s)."""
    tot, k = mp.mpf(0), 1
    while True:
        t = mp.log(R(s * 2 ** (k - 1))) / 2 ** k
        tot += t
        if abs(t) < mp.mpf(10) ** (-mp.mp.dps):
            return tot
        k += 1


P3_2 = mp.e ** logP3(2)
Sigma3 = -mp.diff(logP3, 2)                      # sum_{p=3(4)} log p/(p^2-1)
LpL = mp.euler + 2 * mp.log(2) + 3 * mp.log(mp.pi) - 4 * mp.log(mp.gamma(mp.mpf(1) / 4))
C = mp.sqrt(P3_2) / (mp.pi * mp.sqrt(2))
kappa = mp.mpf(1) / 2 - mp.euler / 4 + mp.log(2) / 4 + Sigma3 / 2 + LpL / 4

# log H = sum_p [log(1+1/(2p)) + (1/2) log(1-1/p)] = sum_{k>=2} a_k P(k)
logH = mp.nsum(lambda k: ((-1) ** (k + 1) * mp.power(2, -k) - mp.mpf(1) / 2) / k * mp.primezeta(k), [2, mp.inf])
H = mp.e ** logH
C2 = H / (2 * mp.sqrt(mp.pi))
# T = sum_p log p /((2p+1)(p-1)) = -(1/2) sum_{k>=0} c_k P'(k+2), c_k = (2/3)(1-(-1/2)^(k+1))
T = -mp.nsum(lambda k: (mp.mpf(2) / 3) * (1 - mp.power(-0.5, k + 1)) * mp.diff(mp.primezeta, k + 2), [0, mp.inf]) / 2
kappa2 = mp.mpf(1) / 4 - mp.euler / 4 - T / 4

for name, v in [("P3(2)", P3_2), ("Sigma3", Sigma3), ("L'/L(1,chi4)", LpL), ("C", C), ("kappa", kappa),
                ("H", H), ("C2", C2), ("T", T), ("kappa2", kappa2)]:
    print(f"{name:14s} {mp.nstr(v, 25)}")
