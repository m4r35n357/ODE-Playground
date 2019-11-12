/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-dbg 16 10 0.1 50000 1 1 0 3 2.2 1 .01
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long n, nsteps;
mpfr_t t, x0, y0, z0, a, b, g, m, h, _, *wa, *wb, *w1, *x, *y, *z;

int main (int argc, char **argv) {
    assert(argc == 12);
    // initialize from command arguments
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 6, &y0);
    t_arg(argv, 7, &z0);
    t_arg(argv, 8, &a);
    t_arg(argv, 9, &b);
    t_arg(argv, 10, &g);
    t_arg(argv, 11, &m);

    // initialize the derivative and temporary jets
    x = t_jet(n + 1);
    y = t_jet(n + 1);
    z = t_jet(n + 1);
    wa = t_jet(n);
    wb = t_jet(n);
    mpfr_init_set_ui(_, 1, RND);
    w1 = t_jet_c(n, _);

    // main loop
    t_xyz_output(x0, y0, z0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        mpfr_set(y[0], y0, RND);
        mpfr_set(z[0], z0, RND);
        for (int k = 0; k < n; k++) {
            //  x' = Ax(1 - y) - Bz
            mpfr_sub(wa[k], w1[k], y[k], RND);
            mpfr_fmms(_, a, *t_prod(&_, x, wa, k), b, z[k], RND);
            mpfr_div_ui(x[k + 1], _, k + 1, RND);
            //  y' = - Gy(1 - x^2)
            mpfr_sub(wb[k], w1[k], *t_sqr(&_, x, k), RND);
            mpfr_mul(_, *t_prod(&_, y, wb, k), g, RND);
            mpfr_div_si(y[k + 1], _, - (k + 1), RND);
            //  z' = Mx
            mpfr_mul(_, m, x[k], RND);
            mpfr_div_ui(z[k + 1], _, k + 1, RND);
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
