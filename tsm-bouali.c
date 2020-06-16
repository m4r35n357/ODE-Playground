/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-dbg 9 32 10 0.01 50000 1 1 0 3 2.2 1 .01
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t a, b, g, m, h, _;

    // initialize from command arguments
    assert(argc == 13);
    mpfr_inits(a, b, g, m, _, NULL);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_series(n + 1), y = t_series(n + 1), z = t_series(n + 1);
    t_args(argv, argc, x.jet, y.jet, z.jet, &a, &b, &g, &m);
    series gx2 = t_series(n);

    t_output(x.jet[0], y.jet[0], z.jet[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = Ax(1 - y) - Bz
            mpfr_fmms(_, a, x.jet[k], a, *t_prod(x, y, k), RND);
            mpfr_fms(_, b, z.jet[k], _, RND);
            t_next(x, _, k, NEG);
            //  y' = - Gy(1 - x^2)
            mpfr_mul(gx2.jet[k], g, *t_sqr(x, k), RND);
            mpfr_fms(_, g, y.jet[k], *t_prod(y, gx2, k), RND);
            t_next(y, _, k, NEG);
            //  z' = Mx
            mpfr_mul(_, m, x.jet[k], RND);
            t_next(z, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
