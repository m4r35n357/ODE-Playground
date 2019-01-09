/*
 * Lotka-Volterra (Predator-Prey) System
 *
 * Example: ./tsm-lotka-volterra-dbg 16 10 .01 2001 10 10 1 .5 .05 .02
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, a, b, c, d, h, _, wxy, *cx, *cy;

int main (int argc, char **argv) {
    assert(argc == 11);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    mpfr_inits(_, wxy, NULL);
    mpfr_init_set_str(x, argv[5], BASE, RND);
    mpfr_init_set_str(y, argv[6], BASE, RND);
    mpfr_init_set_str(a, argv[7], BASE, RND);
    mpfr_init_set_str(b, argv[8], BASE, RND);
    mpfr_init_set_str(c, argv[9], BASE, RND);
    mpfr_init_set_str(d, argv[10], BASE, RND);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);

    // main loop
    for (long step = 1; step < nsteps + 1; step++) {
        // print a line of output
        t_line_output(t, 3, x, y, x);

        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        for (int k = 0; k < order; k++) {
            t_product(&wxy, cx, cy, k);
            //  x' = Ax - Cxy
            mpfr_fmms(_, a, cx[k], c, wxy, RND);
            mpfr_div_ui(cx[k + 1], _, k + 1, RND);
            //  y' = Dxy - By
            mpfr_fmms(_, d, wxy, b, cy[k], RND);
            mpfr_div_ui(cy[k + 1], _, k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x, cx, order, h);
        t_horner(&y, cy, order, h);
        mpfr_mul_ui(t, h, step, RND);
    }
    return 0;
}
