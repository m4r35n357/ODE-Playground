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
    mpfr_t x0, y0, z0, a, b, h, _;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    series x = t_jet_c(n + 1, x0), y = t_jet_c(n + 1, y0), z = t_jet_c(n + 1, z0);
    series sax = t_jet(n), say = t_jet(n), saz = t_jet(n);
    series cax = t_jet(n), cay = t_jet(n), caz = t_jet(n);
    series ax = t_jet(n), ay = t_jet(n), az = t_jet(n);
    series tx = t_jet(n), ty = t_jet(n), tz = t_jet(n);
    series sx = t_jet(n), sy = t_jet(n), sz = t_jet(n);

    t_output(x.a[0], y.a[0], z.a[0], h, 0, _);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            mpfr_mul(ax.a[k], x.a[k], a, RND);
            mpfr_mul(ay.a[k], y.a[k], a, RND);
            mpfr_mul(az.a[k], z.a[k], a, RND);
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
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step, _);
    }
    return 0;
}
