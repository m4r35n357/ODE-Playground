/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-dbg 9 32 10 0.01 30000 1 0 0 4.75 1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t a, b, h, _;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &h, &nsteps);
    mpfr_inits(a, b, _, NULL);
    series x = t_series(n + 1), y = t_series(n + 1), z = t_series(n + 1);
    t_args(argv, argc, x.jet, y.jet, z.jet, &a, &b);
    series sax = t_series(n), say = t_series(n), saz = t_series(n);
    series cax = t_series(n), cay = t_series(n), caz = t_series(n);
    series ax = t_series(n), ay = t_series(n), az = t_series(n);
    series tx = t_series(n), ty = t_series(n), tz = t_series(n);
    series sx = t_series(n), sy = t_series(n), sz = t_series(n);

    t_output(x.jet[0], y.jet[0], z.jet[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            mpfr_mul(ax.jet[k], x.jet[k], a, RND);
            mpfr_mul(ay.jet[k], y.jet[k], a, RND);
            mpfr_mul(az.jet[k], z.jet[k], a, RND);
            //  x' = sin(Ay) - Btan(x)
            mpfr_fms(_, b, *t_tan_sec2(tx, sx, x, k, TRIG).a, *t_sin_cos(say, cay, ay, k, TRIG).a, RND);
            t_next(x, _, k, NEG);
            //  y' = sin(Az) - Btan(y)
            mpfr_fms(_, b, *t_tan_sec2(ty, sy, y, k, TRIG).a, *t_sin_cos(saz, caz, az, k, TRIG).a, RND);
            t_next(y, _, k, NEG);
            //  z' = sin(Ax) - Btan(z)
            mpfr_fms(_, b, *t_tan_sec2(tz, sz, z, k, TRIG).a, *t_sin_cos(sax, cax, ax, k, TRIG).a, RND);
            t_next(z, _, k, NEG);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
