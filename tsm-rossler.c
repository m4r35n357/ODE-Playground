/*
 * Rossler System
 *
 * Example: ./tsm-rossler-dbg 16 10 0.01 150001 0.0 -6.78 0.02 .2 .2 5.7
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, a, b, c, h, tmp, *w_b, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 11);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    mpfr_inits(tmp, NULL);
    mpfr_init_set_str(x, argv[5], BASE, RND);
    mpfr_init_set_str(y, argv[6], BASE, RND);
    mpfr_init_set_str(z, argv[7], BASE, RND);
    mpfr_init_set_str(a, argv[8], BASE, RND);
    mpfr_init_set_str(b, argv[9], BASE, RND);
    mpfr_init_set_str(c, argv[10], BASE, RND);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    w_b = t_jet_c(order, b);

    // main loop
    for (long step = 1; step < nsteps + 1; step++) {
        // print a line of output
        t_line_output(t, 3, x, y, z);

        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            //  x' = - y - z
            mpfr_add(tmp, cy[k], cz[k], RND);
            mpfr_div_si(cx[k + 1], tmp, - (k + 1), RND);
            //  y' = x + Ay
            mpfr_fma(tmp, a, cy[k], cx[k], RND);
            mpfr_div_ui(cy[k + 1], tmp, k + 1, RND);
            //  z' = B + z(x - C)
            t_product(&tmp, cz, cx, k);
            mpfr_add(tmp, w_b[k], tmp, RND);
            mpfr_fms(tmp, c, cz[k], tmp, RND);
            mpfr_div_si(cz[k + 1], tmp, - (k + 1), RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x, cx, order, h);
        t_horner(&y, cy, order, h);
        t_horner(&z, cz, order, h);
        mpfr_mul_ui(t, h, step, RND);
    }
    return 0;
}
