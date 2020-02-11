/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-dbg 16 10 .01 50000 .05 -.05 .3 .28713 .1
 *          ./tsm-rf-dbg 32 16 .01 50000 .05 -.05 .3 .105 .1 | ./plotPi3d.py
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, gamma, h, d3, _, x2_1, *a, *b, *c, *w1, *alpha, *x, *y, *z;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &_, &gamma);
    alpha = t_jet_c(n, _);

    // initialize the derivative and temporary jets
    x = t_jet(n + 1);
    y = t_jet(n + 1);
    z = t_jet(n + 1);
    a = t_jet(n);
    b = t_jet(n);
    c = t_jet(n);
    mpfr_set_ui(_, 1, RND);
    w1 = t_jet_c(n, _);
    mpfr_init_set_ui(d3, 3, RND);
    mpfr_init(x2_1);

    // main loop
    t_xyz_output(x0, y0, z0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        mpfr_set(y[0], y0, RND);
        mpfr_set(z[0], z0, RND);
        for (int k = 0; k < n; k++) {
            mpfr_sub(x2_1, *t_sqr(&_, x, k), w1[k], RND);
            //  x' = y(z - 1 + x^2) + Gx
            mpfr_add(a[k], z[k], x2_1, RND);
            mpfr_fma(_, gamma, x[k], *t_prod(&_, y, a, k), RND);
            mpfr_div_ui(x[k + 1], _, k + 1, RND);
            //  y' = x(3z + 1 - x^2) + Gy
            mpfr_fms(b[k], d3, z[k], x2_1, RND);
            mpfr_fma(_, gamma, y[k], *t_prod(&_, x, b, k), RND);
            mpfr_div_ui(y[k + 1], _, k + 1, RND);
            //  z' = -2z(A + xy)
            mpfr_add(c[k], *t_prod(&_, x, y, k), alpha[k], RND);
            mpfr_mul_2ui(_, *t_prod(&_, z, c, k), 1, RND);
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
