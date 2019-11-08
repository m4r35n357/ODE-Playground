/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-dbg 16 10 0.01 30001 1 0 0 4.75 1
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long n, nsteps;
mpfr_t t, x0, y0, z0, a, b, h, _, *tx, *s2x, *ty, *s2y, *tz, *s2z, *ax, *ay, *az,
        *sax, *cax, *say, *cay, *saz, *caz, *x, *y, *z;

int main (int argc, char **argv) {
    assert(argc == 10);
    // initialize from command arguments
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 6, &y0);
    t_arg(argv, 7, &z0);
    t_arg(argv, 8, &a);
    t_arg(argv, 9, &b);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    x = t_jet(n + 1);
    y = t_jet(n + 1);
    z = t_jet(n + 1);
    ax = t_jet(n);
    ay = t_jet(n);
    az = t_jet(n);
    sax = t_jet(n);
    say = t_jet(n);
    saz = t_jet(n);
    cax = t_jet(n);
    cay = t_jet(n);
    caz = t_jet(n);
    tx = t_jet(n);
    ty = t_jet(n);
    tz = t_jet(n);
    s2x = t_jet(n);
    s2y = t_jet(n);
    s2z = t_jet(n);

    // main loop
    t_xyz_output(x0, y0, z0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        mpfr_set(y[0], y0, RND);
        mpfr_set(z[0], z0, RND);
        for (int k = 0; k < n; k++) {
            mpfr_mul(ax[k], x[k], a, RND);
            mpfr_mul(ay[k], y[k], a, RND);
            mpfr_mul(az[k], z[k], a, RND);
            //  x' = sin(Ay) - Btan(x)
            mpfr_fms(_, b, *t_tan_sec2(tx, s2x, x, k, &_, TRIG).a, *t_sin_cos(say, cay, ay, k, &_, TRIG).a, RND);
            mpfr_div_si(x[k + 1], _, - (k + 1), RND);
            //  y' = sin(Az) - Btan(y)
            mpfr_fms(_, b, *t_tan_sec2(ty, s2y, y, k, &_, TRIG).a, *t_sin_cos(saz, caz, az, k, &_, TRIG).a, RND);
            mpfr_div_si(y[k + 1], _, - (k + 1), RND);
            //  z' = sin(Ax) - Btan(z)
            mpfr_fms(_, b, *t_tan_sec2(tz, s2z, z, k, &_, TRIG).a, *t_sin_cos(sax, cax, ax, k, &_, TRIG).a, RND);
            mpfr_div_si(z[k + 1], _, - (k + 1), RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x0, x, n, h);
        t_horner(&y0, y, n, h);
        t_horner(&z0, z, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x0, y0, z0, t);
    }
    return 0;
}
