/* Count A(y) = #{ a <= y : a odd, every prime factor of a is 3 mod 4 }.
   The exceptional set is E = { 2a : a counted by A }, so E(x) = A(x/2).
   Usage: exceptions_count Y   (prints A(y) at checkpoints y = c*10^k <= Y) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define GET(b,i) ((b[(i)>>3]>>((i)&7))&1)
#define SET(b,i) (b[(i)>>3]|=(uint8_t)(1u<<((i)&7)))

int main(int argc, char **argv) {
    uint64_t Y = strtoull(argv[1], 0, 10);
    uint64_t M = Y / 2 + 1;               /* odd a = 2i+1, i < M */
    uint8_t *comp = calloc(M / 8 + 1, 1); /* composite flags for odd numbers */
    uint8_t *bad  = calloc(M / 8 + 1, 1); /* has a prime factor 1 mod 4 */
    SET(comp, 0);                         /* 1 is not prime */
    for (uint64_t i = 1; ; i++) {
        uint64_t p = 2 * i + 1;
        if (p * p > Y) break;
        if (GET(comp, i)) continue;
        for (uint64_t j = (p * p) / 2; j < M; j += p) SET(comp, j);
    }
    for (uint64_t i = 1; i < M; i++) {
        uint64_t p = 2 * i + 1;
        if (p > Y) break;
        if (GET(comp, i) || (p & 3) != 1) continue;
        for (uint64_t j = i; j < M; j += p) SET(bad, j);
    }
    uint64_t cps[64]; int nc = 0;
    for (uint64_t t = 10; t <= Y; t *= 10) {
        uint64_t c[3] = {t, 2 * t, 5 * t};
        for (int k = 0; k < 3; k++) if (c[k] <= Y) cps[nc++] = c[k];
    }
    uint64_t cnt = 0; int ci = 0;
    for (uint64_t i = 0; i < M && ci < nc; i++) {
        uint64_t a = 2 * i + 1;
        while (ci < nc && a > cps[ci]) { printf("%llu %llu\n", (unsigned long long)cps[ci], (unsigned long long)cnt); ci++; }
        if (a > Y) break;
        if (!GET(bad, i)) cnt++;
    }
    while (ci < nc) { printf("%llu %llu\n", (unsigned long long)cps[ci], (unsigned long long)cnt); ci++; }
    return 0;
}
