/*
 * Yu-Wang System
 *
 * Example: ./tsm-yu-wang-dbg 9 32 10 .001 50000 1 0 0 10 40 2 2.5
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t a, b, c, d, h, _;

    // initialize from command arguments
    assert(argc == 13);
    t_stepper(argv, &n, &h, &nsteps);
    mpfr_inits(a, b, c, d, _, NULL);
    series x = t_series(n + 1), y = t_series(n + 1), z = t_series(n + 1);
    t_args(argv, argc, x.jet, y.jet, z.jet, &a, &b, &c, &d);
    series xy = t_series(n), e_xy = t_series(n);

    t_output(x.jet[0], y.jet[0], z.jet[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = A(y - x)
            mpfr_fmms(_, a, y.jet[k], a, x.jet[k], RND);
            t_next(x, _, k, POS);
            //  y' = Bx - cxz
            mpfr_fmms(_, b, x.jet[k], c, *t_prod(x, z, k), RND);
            t_next(y, _, k, POS);
            //  z' = e^(xy) - Dz
            mpfr_set(xy.jet[k], *t_prod(x, y, k), RND);
            mpfr_fms(_, d, z.jet[k], *t_exp(e_xy, xy, k), RND);
            t_next(z, _, k, NEG);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
