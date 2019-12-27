/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-dbg 16 10 0.1 30000 1 0 0 .19
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, b, h, _, *sy, *cy, *sz, *cz, *sx, *cx, *x, *y, *z;

    // initialize from command arguments
    assert(argc == 9);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 6, &y0);
    t_arg(argv, 7, &z0);
    t_arg(argv, 8, &b);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    x = t_jet(n + 1);
    y = t_jet(n + 1);
    z = t_jet(n + 1);
    sx = t_jet(n);
    sy = t_jet(n);
    sz = t_jet(n);
    cx = t_jet(n);
    cy = t_jet(n);
    cz = t_jet(n);

    // main loop
    t_xyz_output(x0, y0, z0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        mpfr_set(y[0], y0, RND);
        mpfr_set(z[0], z0, RND);
        for (int k = 0; k < n; k++) {
            //  x' = sin(y) - Bx
            mpfr_fms(_, b, x[k], *t_sin_cos(sy, cy, y, k, &_, TRIG).a, RND);
            mpfr_div_si(x[k + 1], _, - (k + 1), RND);
            //  y' = sin(z) - By
            mpfr_fms(_, b, y[k], *t_sin_cos(sz, cz, z, k, &_, TRIG).a, RND);
            mpfr_div_si(y[k + 1], _, - (k + 1), RND);
            //  z' = sin(x) - Bz
            mpfr_fms(_, b, z[k], *t_sin_cos(sx, cx, x, k, &_, TRIG).a, RND);
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
