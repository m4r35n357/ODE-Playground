/*
 * Sprott-Jafari System
 *
 * Example: ./tsm-sj-dbg 16 10 .01 10000 0 3.9 .7 8.888 4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, a, b, h, _, __, *w_b, *x, *y, *z;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b);
    mpfr_inits(_, __, NULL);

    // initialize the derivative and temporary jets
    x = t_jet(n + 1);
    y = t_jet(n + 1);
    z = t_jet(n + 1);
    w_b = t_jet_c(n, b);

    // main loop
    t_xyz_output(x0, y0, z0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        mpfr_set(y[0], y0, RND);
        mpfr_set(z[0], z0, RND);
        for (int k = 0; k < n; k++) {
            //  x' = y
            mpfr_div_ui(x[k + 1], y[k], k + 1, RND);
            //  y' = yz - x
            mpfr_sub(_, *t_prod(&_, y, z, k), x[k], RND);
            mpfr_div_ui(y[k + 1], _, k + 1, RND);
            //  z' = z - ax^2 - y^2 - b
            mpfr_sub(__, z[k], w_b[k], RND);
            mpfr_sub(__, __, *t_sqr(&_, y, k), RND);
            mpfr_fma(_, *t_sqr(&_, x, k), a, __, RND);
            mpfr_div_ui(z[k + 1], _, k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x0, x, n, h);
        t_horner(&y0, y, n, h);
        t_horner(&z0, z, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x0, y0, z0, t);
    }
    return 0;
}
