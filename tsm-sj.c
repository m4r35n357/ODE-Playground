/*
 * Sprott-Jafari System
 *
 * Example: ./tsm-sj-dbg 16 10 .01 10000 0 3.9 .7 8.888 4
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, a, b, h, _, __, *w_b, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 10);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x);
    t_arg(argv, 6, &y);
    t_arg(argv, 7, &z);
    t_arg(argv, 8, &a);
    t_arg(argv, 9, &b);
    mpfr_inits(_, __, NULL);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    w_b = t_jet_c(order, b);

    // main loop
    t_xyz_output(x, y, z, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            //  x' = y
            mpfr_div_ui(cx[k + 1], cy[k], k + 1, RND);
            //  y' = yz - x
            mpfr_sub(_, *t_prod(&_, cy, cz, k), cx[k], RND);
            mpfr_div_ui(cy[k + 1], _, k + 1, RND);
            //  z' = z - ax^2 - y^2 - b
            mpfr_sub(__, cz[k], w_b[k], RND);
            mpfr_sub(__, __, *t_sqr(&_, cy, k), RND);
            mpfr_fma(_, *t_sqr(&_, cx, k), a, __, RND);
            mpfr_div_ui(cz[k + 1], _, k + 1, RND);
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
