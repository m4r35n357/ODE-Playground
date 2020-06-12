/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-dbg 9 32 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t sigma, rho, beta, h, _;

    // initialize from command arguments
    assert(argc == 13);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_jet(n + 1), y = t_jet(n + 1), z = t_jet(n + 1);
    t_args(argv, argc, x.a, y.a, z.a, &sigma, &rho, &beta, &_);
    mpfr_div(beta, beta, _, RND);

    t_output(x.a[0], y.a[0], z.a[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = S(y - x)
            mpfr_fmms(_, sigma, y.a[k], sigma, x.a[k], RND);
            t_next(x, _, k, POS);
            //  y' = x(R - z) - y
            mpfr_fms(_, x.a[k], rho, *t_prod(x, z, k), RND);
            mpfr_sub(_, _, y.a[k], RND);
            t_next(y, _, k, POS);
            //  z' = xy - Bz
            mpfr_fms(_, beta, z.a[k], *t_prod(x, y, k), RND);
            t_next(z, _, k, NEG);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
