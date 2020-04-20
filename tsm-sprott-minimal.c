/*
 * Sprott Minimal System
 *
 * Example: ./tsm-sprott-minimal-dbg 32 4 0.1 10000 .02 0 0 2.017
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t x0, y0, z0, a, h, _;

    // initialize from command arguments
    assert(argc == 9);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);

    t_output(x[0], y[0], z[0], h, 0, _);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = y
            t_next(x, y[k], k, POS);
            //  y' = z
            t_next(y, z[k], k, POS);
            //  z' = - az + y^2 - x
            mpfr_fma(_, a, z[k], x[k], RND);
            mpfr_sub(_, *t_sqr(y, k), _, RND);
            t_next(z, _, k, POS);
        }

        // sum the series using Horner's method and advance one step
        t_horner(x, n, h);
        t_horner(y, n, h);
        t_horner(z, n, h);
        t_output(x[0], y[0], z[0], h, step, _);
    }
    return 0;
}
