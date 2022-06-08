/*
 * Parameter (L, Q) generation for Kerr spacetime
 *
 * Example:  ./h-kerr-gen-light-dbg 3 .8
 *
 ./h-kerr-gen-light-dbg $(yad --title="Step Analysis" --form --separator=" " --align=right \
    --field="r" \
    --field="spin":NUM \
    -- "3" '0.8!0.0..1.0!0.1!1')
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "real.h"

int main(int argc, char **argv) {
    assert(argc == 3);
    real r = strtold(argv[1], NULL);
    real a = strtold(argv[2], NULL);
    real r3 = r * r * r;
    real r2 = r * r;
    real a2 = a * a;
    real L = - (r3 - 3.0L * r2 + a2 * r + a2) / (a * (r - 1.0L));
    real Q = - r3 * (r3 - 6.0L * r2 + 9.0L * r - 4.0L * a2) / (a2 * (r - 1.0L) * (r - 1.0L));
    real r_min = 2.0L * (1.0L + cosl(2.0L / 3.0L * acosl(- fabsl(a))));
    real r_max = 2.0L * (1.0L + cosl(2.0L / 3.0L * acosl(fabsl(a))));
    fprintf(stdout, "\n");
    fprintf(stdout, "a: %.3Lf  M: 1.0  min R = %.3Lf  max R = %.3Lf\nL = %.18Lf  Q = %.18Lf\n", a, r_min, r_max, L, Q);
    if (r < r_min || r > r_max) {
        fprintf(stdout, "r is out of range!\n");
    } else {
        fprintf(stdout, "\n");
        fprintf(stdout, "Simulate:\n");
        fprintf(stdout, "./h-kerr-dbg 6 8 .01 10000 0 .8 0.0 1.0 %La 1.0 %La %.3Lf 0.0 >/tmp/$USER/data\n", L, Q, r);
        fprintf(stdout, "\n");
        fprintf(stdout, "Generate ICs:\n");
        fprintf(stdout, "./h-kerr-dbg 15 8 .01 0 2 .8 0.0 1.0 %La 1.0 %La %.3Lf 0.0\n", L, Q, r);
        fprintf(stdout, "\n");
    }
}
