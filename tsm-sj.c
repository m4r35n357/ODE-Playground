/*
 * Sprott-Jafari System
 *
 * Example: ./tsm-sj-dbg 16 10 .01 100001 0 3.9 .7 8.888 4
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, a, b, h, tmp, tmp2, *w_b, *cx, *cy, *cz;

int main (int argc, char **argv) {
    // initialize from command arguments
    t_stepper(argc, argv, &order, &t, &h, &nsteps);
    mpfr_inits(tmp, tmp2, NULL);
    mpfr_init_set_str(x, argv[5], BASE, RND);
    mpfr_init_set_str(y, argv[6], BASE, RND);
    mpfr_init_set_str(z, argv[7], BASE, RND);
    mpfr_init_set_str(a, argv[8], BASE, RND);
    mpfr_init_set_str(b, argv[9], BASE, RND);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    w_b = t_jet_constant(order, b);

    // main loop
    for (long step = 1; step < nsteps + 1; step++) {
        // print a line of output
        t_line_output(t, 3, x, y, z);

        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            //  x' = y
            mpfr_div_ui(cx[k + 1], cy[k], k + 1, RND);
            //  y' = yz - x
            t_product(&tmp, cy, cz, k);
            mpfr_sub(tmp, tmp, cx[k], RND);
            mpfr_div_ui(cy[k + 1], tmp, k + 1, RND);
            //  z' = z - ax^2 - y^2 - b
            t_square(&tmp2, cy, k);
            mpfr_sub(tmp, cz[k], w_b[k], RND);
            mpfr_sub(tmp2, tmp, tmp2, RND);
            t_square(&tmp, cx, k);
            mpfr_fma(tmp, tmp, a, tmp2, RND);
            mpfr_div_ui(cz[k + 1], tmp, k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x, cx, order, h);
        t_horner(&y, cy, order, h);
        t_horner(&z, cz, order, h);
        mpfr_mul_ui(t, h, step, RND);
    }
    return 0;
}
