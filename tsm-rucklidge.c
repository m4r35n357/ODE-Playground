/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-dbg 9 32 10 0.01 15000 1 0 0 6.7 2
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t x0, y0, z0, alpha, kappa, h, _;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &alpha, &kappa);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    series x = t_jet_c(n + 1, x0), y = t_jet_c(n + 1, y0), z = t_jet_c(n + 1, z0);

    t_output(x.a[0], y.a[0], z.a[0], h, 0, _);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = ay - kx - yz
            mpfr_fmms(_, alpha, y.a[k], kappa, x.a[k], RND);
            mpfr_sub(_, _, *t_prod(y, z, k), RND);
            t_next(x, _, k, POS);
            //  y' = x
            t_next(y, x.a[k], k, POS);
            //  z' = y^2 - z
            mpfr_sub(_, *t_sqr(y, k), z.a[k], RND);
            t_next(z, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step, _);
    }
    return 0;
}
