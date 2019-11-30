/*
 * Yu-Wang System
 *
 * Example: ./tsm-yu-wang-dbg 16 10 .001 50000 1 0 0 10 40 2 2.5
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, a, b, c, d, h, _, *xy, *e_xy, *x, *y, *z;

    // initialize from command arguments
    assert(argc == 12);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 6, &y0);
    t_arg(argv, 7, &z0);
    t_arg(argv, 8, &a);
    t_arg(argv, 9, &b);
    t_arg(argv, 10, &c);
    t_arg(argv, 11, &d);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    x = t_jet(n + 1);
    y = t_jet(n + 1);
    z = t_jet(n + 1);
    xy = t_jet(n);
    e_xy = t_jet(n);

    // main loop
    t_xyz_output(x0, y0, z0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        mpfr_set(y[0], y0, RND);
        mpfr_set(z[0], z0, RND);
        for (int k = 0; k < n; k++) {
            //  x' = A(y - x)
            mpfr_fmms(_, a, y[k], a, x[k], RND);
            mpfr_div_ui(x[k + 1], _, k + 1, RND);
            //  y' = Bx - cxz
            mpfr_mul(_, c, *t_prod(&_, x, z, k), RND);
            mpfr_fms(_, b, x[k], _, RND);
            mpfr_div_ui(y[k + 1], _, k + 1, RND);
            //  z' = e^(xy) - Dz
            mpfr_set(xy[k], *t_prod(&_, x, y, k), RND);
            mpfr_fms(_, d, z[k], *t_exp(e_xy, xy, k, &_), RND);
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
