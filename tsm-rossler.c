/*
 * Rossler System
 *
 * Example: ./tsm-rossler-dbg 9 32 10 0.01 50000 0.0 -6.78 0.02 .2 .2 5.7
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t a, c, h, _;

    // initialize from command arguments
    assert(argc == 12);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_jet(n + 1), y = t_jet(n + 1), z = t_jet(n + 1), b = t_jet(n);
    t_args(argv, argc, x.a, y.a, z.a, &a, b.a, &c);
    mpfr_init(_);

    t_output(x.a[0], y.a[0], z.a[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = - y - z
            mpfr_add(_, y.a[k], z.a[k], RND);
            t_next(x, _, k, NEG);
            //  y' = x + Ay
            mpfr_fma(_, a, y.a[k], x.a[k], RND);
            t_next(y, _, k, POS);
            //  z' = B + z(x - C)
            mpfr_fms(_, c, z.a[k], *t_prod(z, x, k), RND);
            mpfr_sub(_, b.a[k], _, RND);
            t_next(z, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
