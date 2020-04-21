/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-dbg 32 10 0.01 50000 1 1 0 3 2.2 1 .01
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t x0, y0, z0, a, b, g, m, h, _, *gx2;

    // initialize from command arguments
    assert(argc == 12);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b, &g, &m);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);
    gx2 = t_jet(n);

    t_output(x[0], y[0], z[0], h, 0, _);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = Ax(1 - y) - Bz
            mpfr_fmms(_, a, x[k], a, *t_prod(x, y, k), RND);
            mpfr_fms(_, b, z[k], _, RND);
            t_next(x, _, k, NEG);
            //  y' = - Gy(1 - x^2)
            mpfr_mul(gx2[k], g, *t_sqr(x, k), RND);
            mpfr_fms(_, g, y[k], *t_prod(y, gx2, k), RND);
            t_next(y, _, k, NEG);
            //  z' = Mx
            mpfr_mul(_, m, x[k], RND);
            t_next(z, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, n, h), *t_horner(y, n, h), *t_horner(z, n, h), h, step, _);
    }
    return 0;
}
