/*
 * Rossler System
 *
 * Example: ./tsm-rossler-dbg 16 10 0.01 150001 0.0 -6.78 0.02 .2 .2 5.7
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, a, b, c, h, _, *wb, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 11);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x);
    t_arg(argv, 6, &y);
    t_arg(argv, 7, &z);
    t_arg(argv, 8, &a);
    t_arg(argv, 9, &b);
    t_arg(argv, 10, &c);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    wb = t_jet_c(order, b);

    // main loop
    t_xyz_output(x, y, z, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            //  x' = - y - z
            mpfr_add(_, cy[k], cz[k], RND);
            mpfr_div_si(cx[k + 1], _, - (k + 1), RND);
            //  y' = x + Ay
            mpfr_fma(_, a, cy[k], cx[k], RND);
            mpfr_div_ui(cy[k + 1], _, k + 1, RND);
            //  z' = B + z(x - C)
            t_prod(&_, cz, cx, k);
            mpfr_add(_, wb[k], _, RND);
            mpfr_fms(_, c, cz[k], _, RND);
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
