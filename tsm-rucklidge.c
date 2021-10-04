/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-dbg 9 32 10 0.01 15000 1 0 0 6.7 2
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t alpha, kappa, h, _;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &h, &nsteps);
    mpfr_inits(alpha, kappa, _, NULL);
    series x = t_series(n + 1), y = t_series(n + 1), z = t_series(n + 1);
    t_args(argv, argc, x.jet, y.jet, z.jet, &alpha, &kappa);

    t_output(x.jet[0], y.jet[0], z.jet[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = ay - kx - yz
            mpfr_fmms(_, alpha, y.jet[k], kappa, x.jet[k], RND);
            mpfr_sub(_, _, *t_prod(y, z, k), RND);
            t_next(x, _, k, POS);
            //  y' = x
            t_next(y, x.jet[k], k, POS);
            //  z' = y^2 - z
            mpfr_sub(_, *t_sqr(y, k), z.jet[k], RND);
            t_next(z, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
