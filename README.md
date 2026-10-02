# On the ratio of surjective group homomorphisms to ring homomorphisms between cyclic rings

Priyabrata Mandal (corresponding author, MANIT Bhopal), Deep Bhattacharjee and [Third Author Name].

Manuscript prepared for *AIMS Mathematics*. The compiled paper is
[`paper/main.pdf`](paper/main.pdf); the source is [`paper/main.tex`](paper/main.tex).

## Results

For `n | m` there are `phi(n)` surjective group homomorphisms `Z_m -> Z_n` and
`2^omega(n)` ring homomorphisms. With `rho(n) = phi(n) / 2^omega(n)`:

1. `rho(n) = 2^eps(n) |U_n^2|`, where `U_n^2` is the group of squares in
   `(Z/nZ)^x` and `eps(n)` is `0, -1, 0, 1` according as `2^0, 2^1, 2^2, 2^(>=3)`
   exactly divides `n`. Hence `rho(n)` is an integer unless `n = 2a` with every
   prime factor of `a` congruent to 3 mod 4.
2. For odd `n` the ring homomorphisms correspond canonically to the subfields of
   `Q(zeta_n)` of degree at most 2, and `rho(n) = [Q(zeta_n) : M_n]` with `M_n`
   the maximal multiquadratic subfield.
3. The number of `n <= x` with `rho(n)` not an integer is
   `C x (log x)^(-1/2) (1 + kappa/log x + O((log x)^-2))`, with
   `C = 0.2432599441929...` and `kappa = 0.7047534517059...`.
4. `sum_{n <= x} rho(n) = C2 x^2 (log x)^(-1/2) (1 + kappa2/log x + O((log x)^-2))`,
   with `C2 = 0.2290907914889...` and `kappa2 = 0.0248264240338...`.

## Reproducing the numbers

Requirements: a C compiler, Python 3 with `mpmath`, `numpy` and `matplotlib`, and a
TeX distribution (with `cm-super`) for the figures and the paper.

| Command | What it does | Time |
|---|---|---|
| `python3 code/verify_algebra.py` | Brute-force check of Lemmas 2.1-2.3, Theorem 1.1 (`n <= 10^4`) and Theorem 1.2 (odd `n <= 1000`, Gauss sums for `n <= 120`) | seconds |
| `python3 code/constants.py` | `C`, `kappa`, `C2`, `kappa2` and the auxiliary constants to 25 digits | ~3 min |
| `gcc -O2 -o exc code/exceptions_count.c -lm && ./exc 10000000000` | Exact counts `A(y)` of odd `a <= y` built from primes 3 mod 4; `E(2y) = A(y)` | ~2 min |
| `gcc -O2 -o rho_sum code/rho_sum.c -lm && ./rho_sum 10000000000` | Exact values of `2 * sum_{n<=x} rho(n)` | ~9 min |
| `python3 code/tables.py` | Rows of Tables 1-3 from `data/` | seconds |
| `python3 code/figures.py` | Figures 1-4 into `paper/figures/` | seconds |
| `cd paper && latexmk -pdf main.tex` | The paper | seconds |
| `sh paper/make_arxiv.sh` | arXiv source bundle `paper/arxiv_source.zip` | seconds |

The outputs of the two C programs are stored in [`data/odd_3mod4_counts.txt`](data/odd_3mod4_counts.txt)
(columns `y`, `A(y)`) and [`data/rho_sums.txt`](data/rho_sums.txt) (columns `x`, `2 * sum_{n<=x} rho(n)`).

## Citation

See [`CITATION.cff`](CITATION.cff).
