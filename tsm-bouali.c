/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-dbg 16 10 0.1 50000 1 1 0 3 2.2 1 .01
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, a, b, g, m, h, _, *wa, *wb, *w1;

    // initialize from command arguments
    assert(argc == 12);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b, &g, &m);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);
    wa = t_jet(n);
    wb = t_jet(n);
    mpfr_init_set_ui(_, 1, RND);
    w1 = t_jet_c(n, _);

    // main loop
    t_xyz_output(x[0], y[0], z[0], t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = Ax(1 - y) - Bz
            mpfr_sub(wa[k], w1[k], y[k], RND);
            mpfr_fmms(_, a, *t_prod(x, wa, k), b, z[k], RND);
            mpfr_div_ui(x[k + 1], _, k + 1, RND);
            //  y' = - Gy(1 - x^2)
            mpfr_sub(wb[k], w1[k], *t_sqr(x, k), RND);
            mpfr_mul(_, *t_prod(y, wb, k), g, RND);
            mpfr_div_si(y[k + 1], _, - (k + 1), RND);
            //  z' = Mx
            mpfr_mul(_, m, x[k], RND);
            mpfr_div_ui(z[k + 1], _, k + 1, RND);
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
