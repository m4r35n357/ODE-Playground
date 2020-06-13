/*
 * Genesio-Tesi System
 *
 * Example: ./tsm-genesio-tesi-dbg 9 32 10 0.01 150000 .1 .1 .1 .44 1.1 1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t a, b, c, h, _;

    // initialize from command arguments
    assert(argc == 12);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_series(n + 1), y = t_series(n + 1), z = t_series(n + 1);
    t_args(argv, argc, x.jet, y.jet, z.jet, &a, &b, &c);
    mpfr_init(_);

    t_output(x.jet[0], y.jet[0], z.jet[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = y
            t_next(x, y.jet[k], k, POS);
            //  y' = z
            t_next(y, z.jet[k], k, POS);
            //  z' = x^2 - Cx - By - Az
            mpfr_fms(_, c, x.jet[k], *t_sqr(x, k), RND);
            mpfr_fma(_, b, y.jet[k], _, RND);
            mpfr_fma(_, a, z.jet[k], _, RND);
            t_next(z, _, k, NEG);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
