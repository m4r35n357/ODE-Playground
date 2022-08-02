/*
 * Parameter (L, Q) generation for Kerr spacetime
 *
 * Example:  ./h-kerr-gen-light-dbg 3 .8
 *
 ./h-kerr-gen-light-dbg $(yad --columns=2 --title="Generate Parameters (light)" --form --separator=" " --align=right \
    --field="r":NUM --field="BH spin":NUM \
    -- '3.0!1.0..4.0!0.1!1' '0.8!-1.0..1.0!0.1!1')
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
        fprintf(stdout, "./h-kerr-dbg 6 8 .01 10000 0 %.3Lf 0.0 1.0 %La 1.0 %La %.3Lf 0.0 >/tmp/$USER/data\n", a, L, Q, r);
        fprintf(stdout, "\n");
        fprintf(stderr, "./h-kerr-gl $(yad --columns=2 --title='Kerr Light Orbit GL' --form --separator=' ' --align=right ");
        fprintf(stderr, "--field='Display Mode':CB ");
        fprintf(stderr, "--field='Order':NUM ");
        fprintf(stderr, "--field='Step Size':NUM ");
        fprintf(stderr, "--field='Steps':NUM ");
        fprintf(stderr, "--field='Track Length':NUM ");
        fprintf(stderr, "--field='BH spin':NUM ");
        fprintf(stderr, "--field='particle mass':RO ");
        fprintf(stderr, "--field='particle energy':RO ");
        fprintf(stderr, "--field='particle momentum' ");
        fprintf(stderr, "--field='momentum factor':RO ");
        fprintf(stderr, "--field='Carter constant' ");
        fprintf(stderr, "--field='r0' ");
        fprintf(stderr, "--field='theta0' ");
        fprintf(stderr, "-- '0!1!2' '4!2..10!1' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' '1000!1..100000!1' "),
        fprintf(stderr, "'%.3Lf!-1.0..1.0!0.1!1' 0.0 1.0 %.9Le 1.0 %.9Le %.3Lf 0.0)\n", a, L, Q, r);
        fprintf(stderr, "\n");
        fprintf(stdout, "Generate ICs:\n");
        fprintf(stdout, "./h-kerr-dbg 15 8 .01 0 2 %.3Lf 0.0 1.0 %La 1.0 %La %.3Lf 0.0\n", a, L, Q, r);
        fprintf(stdout, "\n");
    }
}
