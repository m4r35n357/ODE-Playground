/*
 * Yu-Wang System
 *
 * Example: ./tsm-yu-wang-dbg 32 10 .001 50000 1 0 0 10 40 2 2.5
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, a, b, c, d, h, _;

    // initialize from command arguments
    assert(argc == 12);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b, &c, &d);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);
    mpfr_t *xy = t_jet(n);
    mpfr_t *e_xy = t_jet(n);

    // main loop
    t_xyz_output(x[0], y[0], z[0], t);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = A(y - x)
            mpfr_fmms(_, a, y[k], a, x[k], RND);
            t_next(x, _, k, POS);
            //  y' = Bx - cxz
            mpfr_fmms(_, b, x[k], c, *t_prod(x, z, k), RND);
            t_next(y, _, k, POS);
            //  z' = e^(xy) - Dz
            mpfr_set(xy[k], *t_prod(x, y, k), RND);
            mpfr_fms(_, d, z[k], *t_exp(e_xy, xy, k), RND);
            t_next(z, _, k, NEG);
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
