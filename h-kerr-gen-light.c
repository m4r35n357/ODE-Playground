/*
 * Parameter (L, Q) generation for Kerr spacetime
 *
 * Example:  ./h-kerr-gen-light-std 3 .8
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "real.h"

int main(int argc, char **argv) { (void)argc;
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
    fprintf(stdout, "\na: %.3Lf  M: 1.0  min R = %.3Lf  max R = %.3Lf\nL = %.18Lf  Q = %.18Lf\n", a, r_min, r_max, L, Q);
    if (r < r_min || r > r_max) {
        fprintf(stdout, "r is out of range!\n");
    } else {
        fprintf(stdout, "\nSimulate:\n");
        fprintf(stdout, "./h-kerr-std 6 8 .01 10000 0 %.3Lf 0.0 1.0 %La 1.0 %La %.3Lf 0.0 >/tmp/$USER/data\n", a, L, Q, r);
        fprintf(stderr, "\n./h-kerr-gl $(yad --columns=2 --title='Kerr Light Orbit GL' --form --separator=' ' --align=right ");
        fprintf(stderr, "--field='Display Mode':CB --field='Order':NUM --field='Step Size':NUM --field='Steps':NUM ");
        fprintf(stderr, "--field='Track Length':NUM --field='BH spin':NUM --field='particle mass':RO ");
        fprintf(stderr, "--field='energy':RO --field='momentum' --field='momentum factor':RO --field='Carter constant' ");
        fprintf(stderr, "--field='r0' --field='theta0' ");
        fprintf(stderr, "-- '0!1!2' '4!2..10!2' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' '1000!1..100000!1' "),
        fprintf(stderr, "'%.3Lf!-1.0..1.0!0.1!1' 0.0 1.0 %.9Le 1.0 %.9Le %.3Lf 0.0)\n", a, L, Q, r);
        fprintf(stdout, "\nGenerate ICs:\n");
        fprintf(stdout, "./h-kerr-std 15 8 .01 0 2 %.3Lf 0.0 1.0 %La 1.0 %La %.3Lf 0.0\n\n", a, L, Q, r);
    }
}
