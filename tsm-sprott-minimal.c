/*
 * Sprott Minimal System
 *
 * Example: ./tsm-sprott-minimal-dbg 9 32 4 0.1 10000 .02 0 0 2.017
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
    assert(argc == 10);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    series x = t_jet_c(n + 1, x0), y = t_jet_c(n + 1, y0), z = t_jet_c(n + 1, z0);

    t_output(x.a[0], y.a[0], z.a[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = y
            t_next(x, y.a[k], k, POS);
            //  y' = z
            t_next(y, z.a[k], k, POS);
            //  z' = - az + y^2 - x
            mpfr_fma(_, a, z.a[k], x.a[k], RND);
            mpfr_sub(_, *t_sqr(y, k), _, RND);
            t_next(z, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
