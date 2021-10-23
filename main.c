/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long dp = strtol(argv[1], NULL, BASE); assert(dp >= 1 && dp <= 99);
    long n = strtol(argv[2], NULL, BASE); assert(n >= 2 && n <= 64);
    real h = strtold(argv[3], NULL); assert(h > 0.0L);
    long steps = strtol(argv[4], NULL, BASE); assert(steps >= 1 && steps <= 1000000);

    tsm(argc, argv, dp, n, h, steps, strtold(argv[5], NULL), strtold(argv[6], NULL), strtold(argv[7], NULL));

    return 0;
}
