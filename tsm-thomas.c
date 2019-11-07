/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-dbg 16 10 0.1 30001 1 0 0 .19
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x0, y0, z0, b, h, _, *sy, *cy, *sz, *cz, *sx, *cx, *x, *y, *z;

int main (int argc, char **argv) {
    assert(argc == 9);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 6, &y0);
    t_arg(argv, 7, &z0);
    t_arg(argv, 8, &b);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    x = t_jet(order + 1);
    y = t_jet(order + 1);
    z = t_jet(order + 1);
    sx = t_jet(order);
    sy = t_jet(order);
    sz = t_jet(order);
    cx = t_jet(order);
    cy = t_jet(order);
    cz = t_jet(order);

    // main loop
    t_xyz_output(x0, y0, z0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        mpfr_set(y[0], y0, RND);
        mpfr_set(z[0], z0, RND);
        for (int k = 0; k < order; k++) {
            //  x' = sin(y) - Bx
            mpfr_fms(_, x[k], b, t_sin_cos(sy, cy, y, k, &_, TRIG).a[k], RND);
            mpfr_div_si(x[k + 1], _, - (k + 1), RND);
            //  y' = sin(z) - By
            mpfr_fms(_, y[k], b, t_sin_cos(sz, cz, z, k, &_, TRIG).a[k], RND);
            mpfr_div_si(y[k + 1], _, - (k + 1), RND);
            //  z' = sin(x) - Bz
            mpfr_fms(_, z[k], b, t_sin_cos(sx, cx, x, k, &_, TRIG).a[k], RND);
            mpfr_div_si(z[k + 1], _, - (k + 1), RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x0, x, order, h);
        t_horner(&y0, y, order, h);
        t_horner(&z0, z, order, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x0, y0, z0, t);
    }
    return 0;
}
