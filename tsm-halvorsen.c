/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-dbg 9 32 10 .01 10000 1 0 0 1.4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t x0, y0, z0, a, h, _, w4x, w4y, w4z;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a);
    mpfr_inits(_, w4x, w4y, w4z, NULL);

    // initialize the derivative and temporary jets
    series x = t_jet_c(n + 1, x0), y = t_jet_c(n + 1, y0), z = t_jet_c(n + 1, z0);

    t_output(x.a[0], y.a[0], z.a[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            mpfr_mul_2ui(w4x, x.a[k], 2, RND);
            mpfr_mul_2ui(w4y, y.a[k], 2, RND);
            mpfr_mul_2ui(w4z, z.a[k], 2, RND);
            //  x' = - Ax - 4y - 4z - y^2
            mpfr_fma(_, a, x.a[k], *t_sqr(y, k), RND);
            mpfr_add(_, w4y, _, RND);
            mpfr_add(_, w4z, _, RND);
            t_next(x, _, k, NEG);
            //  y' = - Ay - 4z - 4x - z^2
            mpfr_fma(_, a, y.a[k], *t_sqr(z, k), RND);
            mpfr_add(_, w4z, _, RND);
            mpfr_add(_, w4x, _, RND);
            t_next(y, _, k, NEG);
            //  z' = - Az - 4x - 4y - x^2
            mpfr_fma(_, a, z.a[k], *t_sqr(x, k), RND);
            mpfr_add(_, w4x, _, RND);
            mpfr_add(_, w4y, _, RND);
            t_next(z, _, k, NEG);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
