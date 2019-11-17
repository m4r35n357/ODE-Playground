/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-dbg 16 10 .01 100000 1 0 0 1.4
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, a, h, _, w4x, w4y, w4z, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 9);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x);
    t_arg(argv, 6, &y);
    t_arg(argv, 7, &z);
    t_arg(argv, 8, &a);
    mpfr_inits(_, w4x, w4y, w4z, NULL);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);

    // main loop
    t_xyz_output(x, y, z, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            mpfr_mul_2ui(w4x, cx[k], 2, RND);
            mpfr_mul_2ui(w4y, cy[k], 2, RND);
            mpfr_mul_2ui(w4z, cz[k], 2, RND);
            //  x' = - Ax - 4y - 4z - y^2
            mpfr_fma(_, a, cx[k], *t_sqr(&_, cy, k), RND);
            mpfr_add(_, w4y, _, RND);
            mpfr_add(_, w4z, _, RND);
            mpfr_div_si(cx[k + 1], _, - (k + 1), RND);
            //  y' = - Ay - 4z - 4x - z^2
            mpfr_fma(_, a, cy[k], *t_sqr(&_, cz, k), RND);
            mpfr_add(_, w4z, _, RND);
            mpfr_add(_, w4x, _, RND);
            mpfr_div_si(cy[k + 1], _, - (k + 1), RND);
            //  z' = - Az - 4x - 4y - x^2
            mpfr_fma(_, a, cz[k], *t_sqr(&_, cx, k), RND);
            mpfr_add(_, w4x, _, RND);
            mpfr_add(_, w4y, _, RND);
            mpfr_div_si(cz[k + 1], _, - (k + 1), RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x, cx, order, h);
        t_horner(&y, cy, order, h);
        t_horner(&z, cz, order, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x, y, z, t);
    }
    return 0;
}
