/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-dbg 32 10 0.01 30000 1 0 0 4.75 1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, a, b, h, _;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);
    mpfr_t *ax = t_jet(n);
    mpfr_t *ay = t_jet(n);
    mpfr_t *az = t_jet(n);
    mpfr_t *sax = t_jet(n);
    mpfr_t *say = t_jet(n);
    mpfr_t *saz = t_jet(n);
    mpfr_t *cax = t_jet(n);
    mpfr_t *cay = t_jet(n);
    mpfr_t *caz = t_jet(n);
    mpfr_t *tx = t_jet(n);
    mpfr_t *ty = t_jet(n);
    mpfr_t *tz = t_jet(n);
    mpfr_t *s2x = t_jet(n);
    mpfr_t *s2y = t_jet(n);
    mpfr_t *s2z = t_jet(n);

    // main loop
    t_xyz_output(x[0], y[0], z[0], t);
    for (long step = 1; step <= nsteps; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            mpfr_mul(ax[k], x[k], a, RND);
            mpfr_mul(ay[k], y[k], a, RND);
            mpfr_mul(az[k], z[k], a, RND);
            //  x' = sin(Ay) - Btan(x)
            mpfr_fms(_, b, *t_tan_sec2(tx, s2x, x, k, TRIG).a, *t_sin_cos(say, cay, ay, k, TRIG).a, RND);
            mpfr_div_si(x[k + 1], _, - (k + 1), RND);
            //  y' = sin(Az) - Btan(y)
            mpfr_fms(_, b, *t_tan_sec2(ty, s2y, y, k, TRIG).a, *t_sin_cos(saz, caz, az, k, TRIG).a, RND);
            mpfr_div_si(y[k + 1], _, - (k + 1), RND);
            //  z' = sin(Ax) - Btan(z)
            mpfr_fms(_, b, *t_tan_sec2(tz, s2z, z, k, TRIG).a, *t_sin_cos(sax, cax, ax, k, TRIG).a, RND);
            mpfr_div_si(z[k + 1], _, - (k + 1), RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(x, n, h);
        t_horner(y, n, h);
        t_horner(z, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x[0], y[0], z[0], t);
    }
    return 0;
}
