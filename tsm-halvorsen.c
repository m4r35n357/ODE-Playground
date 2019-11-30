/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-dbg 16 10 .01 10000 1 0 0 1.4
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, a, h, _, w4x, w4y, w4z, *x, *y, *z;

    // initialize from command arguments
    assert(argc == 9);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 6, &y0);
    t_arg(argv, 7, &z0);
    t_arg(argv, 8, &a);
    mpfr_inits(_, w4x, w4y, w4z, NULL);

    // initialize the derivative and temporary jets
    x = t_jet(n + 1);
    y = t_jet(n + 1);
    z = t_jet(n + 1);

    // main loop
    t_xyz_output(x0, y0, z0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        mpfr_set(y[0], y0, RND);
        mpfr_set(z[0], z0, RND);
        for (int k = 0; k < n; k++) {
            mpfr_mul_2ui(w4x, x[k], 2, RND);
            mpfr_mul_2ui(w4y, y[k], 2, RND);
            mpfr_mul_2ui(w4z, z[k], 2, RND);
            //  x' = - Ax - 4y - 4z - y^2
            mpfr_fma(_, a, x[k], *t_sqr(&_, y, k), RND);
            mpfr_add(_, w4y, _, RND);
            mpfr_add(_, w4z, _, RND);
            mpfr_div_si(x[k + 1], _, - (k + 1), RND);
            //  y' = - Ay - 4z - 4x - z^2
            mpfr_fma(_, a, y[k], *t_sqr(&_, z, k), RND);
            mpfr_add(_, w4z, _, RND);
            mpfr_add(_, w4x, _, RND);
            mpfr_div_si(y[k + 1], _, - (k + 1), RND);
            //  z' = - Az - 4x - 4y - x^2
            mpfr_fma(_, a, z[k], *t_sqr(&_, x, k), RND);
            mpfr_add(_, w4x, _, RND);
            mpfr_add(_, w4y, _, RND);
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
