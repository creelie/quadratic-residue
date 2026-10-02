/* Exact partial sums S(x) = sum_{n<=x} rho(n), rho(n) = phi(n)/2^omega(n),
   by a segmented sieve. 2*rho(n) is always an integer, so we accumulate
   T(x) = sum 2*rho(n) exactly in 128-bit arithmetic and print T(x) at
   checkpoints x = c*10^k <= X (c = 1, 2, 5).  Usage: rho_sum X */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef unsigned __int128 u128;

static void print_u128(u128 v) {
    char buf[64]; int i = 63; buf[i] = 0;
    if (v == 0) buf[--i] = '0';
    while (v) { buf[--i] = '0' + (int)(v % 10); v /= 10; }
    printf("%s", buf + i);
}

int main(int argc, char **argv) {
    uint64_t X = strtoull(argv[1], 0, 10);
    uint64_t R = (uint64_t)sqrtl((long double)X) + 2;
    char *isc = calloc(R + 1, 1);
    uint64_t *pr = malloc(sizeof(uint64_t) * (R + 1)); uint64_t np = 0;
    for (uint64_t i = 2; i <= R; i++) if (!isc[i]) { pr[np++] = i; for (uint64_t j = i * i; j <= R; j += i) isc[j] = 1; }
    const uint64_t B = 1 << 22;
    uint64_t *rem = malloc(sizeof(uint64_t) * B);
    uint64_t *val = malloc(sizeof(uint64_t) * B); /* 2^omega(n) * (accumulated phi part) tracked as phi and count */
    uint8_t *om = malloc(B);
    uint64_t cps[64]; int nc = 0;
    for (uint64_t t = 10; t <= X; t *= 10) { uint64_t c[3] = {t, 2*t, 5*t}; for (int k = 0; k < 3; k++) if (c[k] <= X) cps[nc++] = c[k]; }
    int ci = 0; u128 T = 0;
    for (uint64_t L = 1; L <= X; L += B) {
        uint64_t H = L + B - 1; if (H > X) H = X;
        uint64_t len = H - L + 1;
        for (uint64_t k = 0; k < len; k++) { rem[k] = L + k; val[k] = 1; om[k] = 0; }
        for (uint64_t q = 0; q < np; q++) {
            uint64_t p = pr[q]; if (p * p > H) break;
            uint64_t s = ((L + p - 1) / p) * p;
            for (uint64_t n = s; n <= H; n += p) {
                uint64_t k = n - L; uint64_t pk = 1;
                while (rem[k] % p == 0) { rem[k] /= p; pk *= p; }
                val[k] *= (pk / p) * (p - 1); om[k]++;
            }
        }
        for (uint64_t k = 0; k < len; k++) {
            if (rem[k] > 1) { val[k] *= rem[k] - 1; om[k]++; }
            /* val = phi(n); 2*rho = phi / 2^(omega-1) */
            uint64_t twice;
            if (om[k] == 0) twice = 2;            /* n = 1: rho(1) = 1 */
            else twice = val[k] >> (om[k] - 1);
            T += twice;
            uint64_t n = L + k;
            while (ci < nc && n == cps[ci]) { printf("%llu ", (unsigned long long)n); print_u128(T); printf("\n"); fflush(stdout); ci++; }
        }
    }
    return 0;
}
