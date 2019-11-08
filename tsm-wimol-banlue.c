/*
 * Wimol-Banlue System
 *
 * Example: ./tsm-wimol-banlue-dbg 16 8 0.1 10001 1 0 0 2.0
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x0, y0, z0, a, b, c, h, _, __, *tx, *s2x, *wa, *x, *y, *z;

int main (int argc, char **argv) {
    assert(argc == 9);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 6, &y0);
    t_arg(argv, 7, &z0);
    t_arg(argv, 8, &a);
    mpfr_inits(_, __, NULL);

    // initialize the derivative and temporary jets
    x = t_jet(order + 1);
    y = t_jet(order + 1);
    z = t_jet(order + 1);
    tx = t_jet(order);
    s2x = t_jet(order);
    wa = t_jet_c(order, a);

    // main loop
    t_xyz_output(x0, y0, z0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        mpfr_set(y[0], y0, RND);
        mpfr_set(z[0], z0, RND);
        for (int k = 0; k < order; k++) {
            //  x' = y - x
            mpfr_sub(_, y[k], x[k], RND);
            mpfr_div_ui(x[k + 1], _, k + 1, RND);
            //  y' = - z * tan(x)
            mpfr_div_si(y[k + 1], *t_prod(&_, z, t_tan_sec2(tx, s2x, x, k, &_, HYP).a, k), - (k + 1), RND);
            //  z' = - A + x * y + |y|
            mpfr_add(_, *t_prod(&_, x, y, k), *t_abs(&__, y, k), RND);
            mpfr_sub(_, _, wa[k], RND);
            mpfr_div_ui(z[k + 1], _, k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x0, x, order, h);
        t_horner(&y0, y, order, h);
        t_horner(&z0, z, order, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x0, y0, z0, t);
    }
    return 0;
}
