"""Figures for the paper. Run from the repository root:
    python3 code/figures.py
Writes paper/figures/fig1_cube.pdf ... fig4_rho_sum.pdf.
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 10,
    "text.latex.preamble": r"\usepackage{amsmath,amssymb}",
})

OUT = os.path.join(os.path.dirname(__file__), "..", "paper", "figures")
DATA = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(OUT, exist_ok=True)

BLUE, ORANGE, GREEN, RED, GREY = "#1f4e79", "#c55a11", "#2e7d32", "#a61b29", "#7f7f7f"

# constants (code/constants.py)
C_E, KAPPA = 0.2432599441929455, 0.7047534517
C_2, KAPPA2 = 0.2290907914889383, 0.0248264240


def legendre(a, p):
    t = pow(a % p, (p - 1) // 2, p)
    return -1 if t == p - 1 else t


# ---------------------------------------------------------------- Figure 1
def fig_cube():
    """Idempotents of Z/105Z on the cube {0,1}^3, labelled by Q(sqrt d_e)."""
    primes = [3, 5, 7]
    star = {3: -3, 5: 5, 7: -7}
    fig = plt.figure(figsize=(5.6, 4.9))
    ax = fig.add_subplot(projection="3d")
    verts = {}
    for b in range(8):
        bits = [(b >> i) & 1 for i in range(3)]
        e = next(x for x in range(105) if all(x % p == bits[i] for i, p in enumerate(primes)))
        d = 1
        for i, p in enumerate(primes):
            if bits[i]:
                d *= star[p]
        verts[tuple(bits)] = (e, d)
    for v in verts:
        for i in range(3):
            if v[i] == 0:
                w = list(v); w[i] = 1; w = tuple(w)
                ax.plot(*zip(v, w), color=[BLUE, ORANGE, GREEN][i], lw=1.6, alpha=0.85)
    for v, (e, d) in verts.items():
        k = sum(v)
        ax.scatter(*v, s=70, color=RED if k == 0 else BLUE, depthshade=False, zorder=5)
        field = r"\mathbb{Q}" if d == 1 else r"\mathbb{Q}(\sqrt{%d})" % d
        off = np.array(v, float) * 0.22 - 0.11
        ax.text(*(np.array(v) + off), r"$e=%d$" "\n" r"$%s$" % (e, field),
                fontsize=8.5, ha="center", va="center")
    ax.set_xlabel(r"$e \bmod 3$", labelpad=-6)
    ax.set_ylabel(r"$e \bmod 5$", labelpad=-6)
    ax.set_zlabel(r"$e \bmod 7$", labelpad=-6)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_ticks([0, 1])
        axis.set_tick_params(pad=-3, labelsize=8)
    ax.set_xlim(-0.25, 1.25); ax.set_ylim(-0.25, 1.25); ax.set_zlim(-0.25, 1.25)
    ax.view_init(elev=22, azim=-58)
    ax.set_box_aspect((1, 1, 1))
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_cube.pdf"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 2
def fig_torus():
    """(Z/35Z)^x = C4 x C6 drawn on a torus, coloured by (u/5, u/7).
    Meridians: fixed u mod 5 = 2^i; parallels: fixed u mod 7 = 3^j."""
    R, r = 3.0, 1.1
    elev, azim = 48, -60
    fig = plt.figure(figsize=(6.2, 4.6))
    ax = fig.add_subplot(projection="3d")
    th, ph = np.meshgrid(np.linspace(0, 2 * np.pi, 90), np.linspace(0, 2 * np.pi, 45))
    X = (R + r * np.cos(ph)) * np.cos(th)
    Y = (R + r * np.cos(ph)) * np.sin(th)
    Z = r * np.sin(ph)
    ax.plot_surface(X, Y, Z, color="#e3eaf2", alpha=0.22, linewidth=0, shade=True)
    s = np.linspace(0, 2 * np.pi, 200)
    for i in range(4):                       # meridian circles
        t = 2 * np.pi * i / 4
        ax.plot((R + r * np.cos(s)) * np.cos(t), (R + r * np.cos(s)) * np.sin(t), r * np.sin(s),
                color=GREY, lw=0.8)
    for j in range(6):                       # parallel circles
        p = 2 * np.pi * j / 6
        ax.plot((R + r * np.cos(p)) * np.cos(s), (R + r * np.cos(p)) * np.sin(s), r * np.sin(p) + 0 * s,
                color=GREY, lw=0.8)
    cols = {(1, 1): BLUE, (1, -1): GREEN, (-1, 1): ORANGE, (-1, -1): RED}
    view = np.array([np.cos(np.radians(elev)) * np.cos(np.radians(azim)),
                     np.cos(np.radians(elev)) * np.sin(np.radians(azim)),
                     np.sin(np.radians(elev))])
    for i in range(4):          # u = 2^i (mod 5)
        for j in range(6):      # u = 3^j (mod 7)
            u = next(x for x in range(1, 35) if x % 5 == pow(2, i, 5) and x % 7 == pow(3, j, 7))
            key = (legendre(u, 5), legendre(u, 7))
            t, p = 2 * np.pi * i / 4, 2 * np.pi * j / 6
            pos = np.array([(R + r * np.cos(p)) * np.cos(t), (R + r * np.cos(p)) * np.sin(t), r * np.sin(p)])
            n = np.array([np.cos(p) * np.cos(t), np.cos(p) * np.sin(t), np.sin(p)])
            front = float(n @ view) > -0.15
            a = 1.0 if front else 0.35
            ax.scatter(*pos, s=58 if front else 34, color=cols[key], alpha=a, depthshade=False,
                       edgecolors="k", linewidths=0.4)
            q = pos + 0.42 * n
            ax.text(*q, str(u), fontsize=8 if front else 6.5, alpha=a, ha="center", va="center")
    labels = {
        (1, 1): r"$(+,+)$: the squares $U_{35}^{2}$",
        (1, -1): r"$(+,-)$",
        (-1, 1): r"$(-,+)$",
        (-1, -1): r"$(-,-)$",
    }
    for k, c in cols.items():
        ax.scatter([], [], [], color=c, edgecolors="k", linewidths=0.4, label=labels[k])
    ax.legend(loc="upper left", fontsize=8, frameon=False,
              title=r"$\bigl(\tfrac{u}{5},\tfrac{u}{7}\bigr)$", title_fontsize=8.5)
    ax.set_axis_off()
    ax.set_box_aspect((1, 1, 0.38))
    ax.view_init(elev=elev, azim=azim)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_torus.pdf"), bbox_inches="tight")
    plt.close(fig)


def read_pairs(name):
    rows = []
    with open(os.path.join(DATA, name)) as f:
        for line in f:
            a, b = line.split()
            rows.append((int(a), int(b)))
    return rows


# ---------------------------------------------------------------- Figure 3
def fig_exceptions():
    rows = read_pairs("odd_3mod4_counts.txt")   # (y, A(y)), E(2y) = A(y)
    xs = np.array([2 * y for y, a in rows if 2 * y >= 1000], float)
    Es = np.array([a for y, a in rows if 2 * y >= 1000], float)
    L = np.log(xs)
    fig, ax = plt.subplots(figsize=(5.6, 3.3))
    ax.plot(np.log10(xs), Es * np.sqrt(L) / xs, "o", ms=4, color=BLUE,
            label=r"$E(x)\sqrt{\log x}/x$ (exact count)")
    t = np.linspace(3, 10.5, 300); Lt = t * math.log(10)
    ax.plot(t, C_E * (1 + KAPPA / Lt), "-", color=ORANGE, lw=1.4,
            label=r"$C(1+\kappa/\log x)$")
    ax.axhline(C_E, color=GREY, ls="--", lw=1, label=r"$C=0.24325\ldots$")
    ax.set_xlabel(r"$\log_{10} x$")
    ax.legend(frameon=False, fontsize=8.5)
    ax.set_xlim(2.8, 10.6)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig3_exceptions.pdf"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Figure 4
def fig_rho_sum():
    rows = read_pairs("rho_sums.txt")           # (x, 2*sum rho)
    xs = np.array([x for x, t in rows if x >= 1000], float)
    Ss = np.array([t / 2 for x, t in rows if x >= 1000], float)
    L = np.log(xs)
    fig, ax = plt.subplots(figsize=(5.6, 3.3))
    ax.plot(np.log10(xs), Ss * np.sqrt(L) / xs ** 2, "o", ms=4, color=BLUE,
            label=r"$\sqrt{\log x}\,x^{-2}\sum_{n\le x}\rho(n)$ (exact sum)")
    t = np.linspace(3, 10.5, 300); Lt = t * math.log(10)
    ax.plot(t, C_2 * (1 + KAPPA2 / Lt), "-", color=ORANGE, lw=1.4,
            label=r"$C_2(1+\kappa_2/\log x)$")
    ax.axhline(C_2, color=GREY, ls="--", lw=1, label=r"$C_2=0.22909\ldots$")
    ax.set_xlabel(r"$\log_{10} x$")
    ax.legend(frameon=False, fontsize=8.5)
    ax.set_xlim(2.8, 10.6)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig4_rho_sum.pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_cube()
    fig_torus()
    fig_exceptions()
    fig_rho_sum()
