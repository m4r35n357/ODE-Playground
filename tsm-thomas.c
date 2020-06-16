/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-dbg 9 32 10 0.1 30000 1 0 0 .19
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t b, h, _;

    // initialize from command arguments
    assert(argc == 10);
    mpfr_inits(b, _, NULL);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_series(n + 1), y = t_series(n + 1), z = t_series(n + 1);
    t_args(argv, argc, x.jet, y.jet, z.jet, &b);
    series sx = t_series(n), sy = t_series(n), sz = t_series(n), cx = t_series(n), cy = t_series(n), cz = t_series(n);

    t_output(x.jet[0], y.jet[0], z.jet[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = sin(y) - Bx
            mpfr_fms(_, b, x.jet[k], *t_sin_cos(sy, cy, y, k, TRIG).a, RND);
            t_next(x, _, k, NEG);
            //  y' = sin(z) - By
            mpfr_fms(_, b, y.jet[k], *t_sin_cos(sz, cz, z, k, TRIG).a, RND);
            t_next(y, _, k, NEG);
            //  z' = sin(x) - Bz
            mpfr_fms(_, b, z.jet[k], *t_sin_cos(sx, cx, x, k, TRIG).a, RND);
            t_next(z, _, k, NEG);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
