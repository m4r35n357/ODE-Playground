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
    mpfr_t x0, y0, z0, a, b, h, _, __;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b);
    mpfr_inits(_, __, NULL);

    // initialize the derivative and temporary jets
    series x = t_jet_c(n + 1, x0), y = t_jet_c(n + 1, y0), z = t_jet_c(n + 1, z0);
    series w_b = t_jet_c(n, b);

    t_output(x.a[0], y.a[0], z.a[0], h, 0, _);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = y
            t_next(x, y.a[k], k, POS);
            //  y' = yz - x
            mpfr_sub(_, *t_prod(y, z, k), x.a[k], RND);
            t_next(y, _, k, POS);
            //  z' = z - ax^2 - y^2 - b
            mpfr_sub(__, z.a[k], w_b.a[k], RND);
            mpfr_sub(__, __, *t_sqr(y, k), RND);
            mpfr_fma(_, *t_sqr(x, k), a, __, RND);
            t_next(z, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step, _);
    }
    return 0;
}
