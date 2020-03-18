/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-dbg 16 10 .01 10000 1 0 0 1.4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, a, h, _, w4x, w4y, w4z;

    // initialize from command arguments
    assert(argc == 9);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a);
    mpfr_inits(_, w4x, w4y, w4z, NULL);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);

    // main loop
    t_xyz_output(x[0], y[0], z[0], t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            mpfr_mul_2ui(w4x, x[k], 2, RND);
            mpfr_mul_2ui(w4y, y[k], 2, RND);
            mpfr_mul_2ui(w4z, z[k], 2, RND);
            //  x' = - Ax - 4y - 4z - y^2
            mpfr_fma(_, a, x[k], *t_sqr(y, k), RND);
            mpfr_add(_, w4y, _, RND);
            mpfr_add(_, w4z, _, RND);
            mpfr_div_si(x[k + 1], _, - (k + 1), RND);
            //  y' = - Ay - 4z - 4x - z^2
            mpfr_fma(_, a, y[k], *t_sqr(z, k), RND);
            mpfr_add(_, w4z, _, RND);
            mpfr_add(_, w4x, _, RND);
            mpfr_div_si(y[k + 1], _, - (k + 1), RND);
            //  z' = - Az - 4x - 4y - x^2
            mpfr_fma(_, a, z[k], *t_sqr(x, k), RND);
            mpfr_add(_, w4x, _, RND);
            mpfr_add(_, w4y, _, RND);
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
