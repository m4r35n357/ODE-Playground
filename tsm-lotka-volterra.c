/*
 * Lotka-Volterra (Predator-Prey) System
 *
 * Example: ./tsm-lotka-volterra-dbg 16 10 .01 2001 10 10 0 1 .5 .05 .02
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, a, b, c, d, h, _, wxy, *cx, *cy;

int main (int argc, char **argv) {
    assert(argc == 12);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x);
    t_arg(argv, 6, &y);
    t_arg(argv, 8, &a);
    t_arg(argv, 9, &b);
    t_arg(argv, 10, &c);
    t_arg(argv, 11, &d);
    mpfr_inits(_, wxy, NULL);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);

    // main loop
    t_xyz_output(x, y, x, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        for (int k = 0; k < order; k++) {
            t_prod(&wxy, cx, cy, k);
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
        t_xyz_output(x, y, x, t);
    }
    return 0;
}