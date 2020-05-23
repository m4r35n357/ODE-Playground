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
    mpfr_t x0, y0, z0, b, h, _;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &b);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    series x = t_jet_c(n + 1, x0), y = t_jet_c(n + 1, y0), z = t_jet_c(n + 1, z0);
    series sx = t_jet(n), sy = t_jet(n), sz = t_jet(n);
    series cx = t_jet(n), cy = t_jet(n), cz = t_jet(n);

    t_output(x.a[0], y.a[0], z.a[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = sin(y) - Bx
            mpfr_fms(_, b, x.a[k], *t_sin_cos(sy, cy, y, k, TRIG).a, RND);
            t_next(x, _, k, NEG);
            //  y' = sin(z) - By
            mpfr_fms(_, b, y.a[k], *t_sin_cos(sz, cz, z, k, TRIG).a, RND);
            t_next(y, _, k, NEG);
            //  z' = sin(x) - Bz
            mpfr_fms(_, b, z.a[k], *t_sin_cos(sx, cx, x, k, TRIG).a, RND);
            t_next(z, _, k, NEG);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
