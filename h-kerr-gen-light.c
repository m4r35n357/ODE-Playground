/*
 * Parameter (L, Q) generation for Kerr spacetime
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "real.h"

int main(int argc, char **argv) { (void)argc;
    PRINT_ARGS(argc, argv);
    CHECK(argc == 3);
    real r = strtold(argv[1], NULL);
    real a = strtold(argv[2], NULL); CHECK(a >= -1.0L && a <= 1.0L);
    real r_min = 2.0L * (1.0L + cosl(2.0L / 3.0L * acosl(- fabsl(a)))); CHECK(r > r_min);
    real r_max = 2.0L * (1.0L + cosl(2.0L / 3.0L * acosl(fabsl(a))));   CHECK(r < r_max);
    real r2 = SQR(r);
    real r3 = r2 * r;
    real a2 = SQR(a);
    real L = - (r3 - 3.0L * r2 + a2 * r + a2) / (a * (r - 1.0L));
    real Q = - r3 * (r3 - 6.0L * r2 + 9.0L * r - 4.0L * a2) / (a2 * (r - 1.0L) * (r - 1.0L));
    fprintf(stdout, "\na: %.3Lf  M: 1.0  min R = %.3Lf  max R = %.3Lf\nL = %.18Lf  Q = %.18Lf\n", a, r_min, r_max, L, Q);
    fprintf(stdout, "\nSimulate:\n");
    fprintf(stdout, "./h-kerr-std 6 8 .01 10000 0 %.3Lf 0.0 1.0 %La 1.0 %La %.3Lf 0.0 >/tmp/$USER/data\n", a, L, Q, r);
    fprintf(stderr, "./h-kerr-gl $(yad --columns=2 --title='Kerr Light Orbit GL' --form --separator=' ' --align=right ");
    fprintf(stderr, "--field='Trail Length':NUM --field='Order':NUM --field='Step Size':NUM --field='Steps':NUM ");
    fprintf(stderr, "--field='BH spin':NUM --field='particle mass':RO ");
    fprintf(stderr, "--field='energy':RO --field='momentum' --field='momentum factor':RO --field='Carter constant' ");
    fprintf(stderr, "--field='r0' --field='theta0' ");
    fprintf(stderr, "-- '2000!1000..10000!1000' '4!2..10!2' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' "),
    fprintf(stderr, "'%.3Lf!-1.0..1.0!0.1!1' 0.0 1.0 %.9Le 1.0 %.9Le %.3Lf 0.0)\n", a, L, Q, r);
}
