
/*
 * Parameter (L, Q) generation for Kerr spacetime
 *
 * Example:  ./h-kerr-gen-light-dbg 3 1 .8
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

typedef long double real;

int main(int argc, char **argv) {
    assert(argc >= 4);
    real r = strtold(argv[1], NULL);
    real m = strtold(argv[2], NULL);
    real a = strtold(argv[3], NULL);
    real r3 = r * r * r;
    real r2 = r * r;
    real a2 = a * a;
    real L = - (r3 - 3.0L * m * r2 + a2 * r + a2 * m) / (a * (r - m));
    real Q = - r3 * (r3 - 6.0L * m * r2 + 9.0L * m * r - 4.0L * a2 * m) / (a2 * (r - m) * (r - m));
    real r_min = 2.0L * m * (1.0L + cosl(2.0L / 3.0L * acosl(- fabsl(a) / m)));
    real r_max = 2.0L * m * (1.0L + cosl(2.0L / 3.0L * acosl(fabsl(a) / m)));
    fprintf(stdout, "%s %.3Lf %s %.3Lf %s %.3Lf %s %.3Lf\n", "a:", a, " M:", m, " min R =", r_min, " max R =", r_max);
    if (r < r_min || r > r_max) {
        fprintf(stdout, "r is out of range!\n");
    } else {
        fprintf(stdout, "./h-kerr-sd-dbg 6 4 .01 10000 .8 1.0 0.0 1.0 %.18Lf 1.0 %.18Lf %.3Lf 0.0 0\n", L, Q, r);
    }
}
