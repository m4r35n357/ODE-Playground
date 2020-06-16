/*
 * Sprott-Jafari System
 *
 * Example: ./tsm-sj-dbg 9 32 10 .01 10000 0 3.9 .7 8.888 4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t a, h, _, __;

    // initialize from command arguments
    assert(argc == 11);
    mpfr_inits(a, _, __, NULL);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_series(n + 1), y = t_series(n + 1), z = t_series(n + 1);
    series w_b = t_series(n);
    t_args(argv, argc, x.jet, y.jet, z.jet, &a, w_b.jet);

    t_output(x.jet[0], y.jet[0], z.jet[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = y
            t_next(x, y.jet[k], k, POS);
            //  y' = yz - x
            mpfr_sub(_, *t_prod(y, z, k), x.jet[k], RND);
            t_next(y, _, k, POS);
            //  z' = z - ax^2 - y^2 - b
            mpfr_sub(__, z.jet[k], w_b.jet[k], RND);
            mpfr_sub(__, __, *t_sqr(y, k), RND);
            mpfr_fma(_, *t_sqr(x, k), a, __, RND);
            t_next(z, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
