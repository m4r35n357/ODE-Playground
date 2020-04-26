/*
 * Genesio-Tesi System
 *
 * Example: ./tsm-genesio-tesi-dbg 32 10 0.01 150000 .1 .1 .1 .44 1.1 1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t x0, y0, z0, a, b, c, h, _;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b, &c);
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
            //  z' = x^2 - Cx - By - Az
            mpfr_fms(_, c, x[k], *t_sqr(x, k), RND);
            mpfr_fma(_, b, y[k], _, RND);
            mpfr_fma(_, a, z[k], _, RND);
            t_next(z, _, k, NEG);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, n, h), *t_horner(y, n, h), *t_horner(z, n, h), h, step, _);
    }
    return 0;
}
