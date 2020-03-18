/*
 * Rossler System
 *
 * Example: ./tsm-rossler-dbg 16 10 0.01 50000 0.0 -6.78 0.02 .2 .2 5.7
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, a, c, h, _;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &_, &c);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);
    mpfr_t *b = t_jet_c(n, _);

    // main loop
    t_xyz_output(x[0], y[0], z[0], t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = - y - z
            mpfr_add(_, y[k], z[k], RND);
            mpfr_div_si(x[k + 1], _, - (k + 1), RND);
            //  y' = x + Ay
            mpfr_fma(_, a, y[k], x[k], RND);
            mpfr_div_ui(y[k + 1], _, k + 1, RND);
            //  z' = B + z(x - C)
            mpfr_add(_, b[k], *t_prod(z, x, k), RND);
            mpfr_fms(_, c, z[k], _, RND);
            mpfr_div_si(z[k + 1], _, - (k + 1), RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&_, x, n, h);
        t_horner(&_, y, n, h);
        t_horner(&_, z, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x[0], y[0], z[0], t);
    }
    return 0;
}
