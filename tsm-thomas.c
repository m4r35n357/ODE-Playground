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
    mpfr_t t, x0, y0, z0, b, h, _;

    // initialize from command arguments
    assert(argc == 9);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &b);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);
    mpfr_t *sx = t_jet(n);
    mpfr_t *sy = t_jet(n);
    mpfr_t *sz = t_jet(n);
    mpfr_t *cx = t_jet(n);
    mpfr_t *cy = t_jet(n);
    mpfr_t *cz = t_jet(n);

    // main loop
    t_xyz_output(x[0], y[0], z[0], t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = sin(y) - Bx
            mpfr_fms(_, b, x[k], *t_sin_cos(sy, cy, y, k, TRIG).a, RND);
            mpfr_div_si(x[k + 1], _, - (k + 1), RND);
            //  y' = sin(z) - By
            mpfr_fms(_, b, y[k], *t_sin_cos(sz, cz, z, k, TRIG).a, RND);
            mpfr_div_si(y[k + 1], _, - (k + 1), RND);
            //  z' = sin(x) - Bz
            mpfr_fms(_, b, z[k], *t_sin_cos(sx, cx, x, k, TRIG).a, RND);
            mpfr_div_si(z[k + 1], _, - (k + 1), RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&_, x, n, h);
        t_horner(&_, y, n, h);
        t_horner(&_, z, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x[0], y[0], z[0], t);
    }
    return 0;
}
