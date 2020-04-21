/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-dbg 32 10 .01 50000 .05 -.05 .3 .28713 .1
 *          ./tsm-rf-dbg 32 16 .01 50000 .05 -.05 .3 .105 .1 | ./plot3d.py
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t x0, y0, z0, gamma, h, d3, _, x2_1;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &_, &gamma);
    mpfr_inits(d3, x2_1, NULL);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);
    mpfr_t *alpha = t_jet_c(n, _);
    mpfr_t *a = t_jet(n);
    mpfr_t *b = t_jet(n);
    mpfr_t *c = t_jet(n);
    mpfr_set_ui(d3, 3, RND);
    mpfr_set_ui(_, 1, RND);
    mpfr_t *w1 = t_jet_c(n, _);

    t_output(x[0], y[0], z[0], h, 0, _);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            mpfr_sub(x2_1, *t_sqr(x, k), w1[k], RND);
            //  x' = y(z - 1 + x^2) + Gx
            mpfr_add(a[k], z[k], x2_1, RND);
            mpfr_fma(_, gamma, x[k], *t_prod(y, a, k), RND);
            t_next(x, _, k, POS);
            //  y' = x(3z + 1 - x^2) + Gy
            mpfr_fms(b[k], d3, z[k], x2_1, RND);
            mpfr_fma(_, gamma, y[k], *t_prod(x, b, k), RND);
            t_next(y, _, k, POS);
            //  z' = -2z(A + xy)
            mpfr_add(c[k], *t_prod(x, y, k), alpha[k], RND);
            mpfr_mul_2ui(_, *t_prod(z, c, k), 1, RND);
            t_next(z, _, k, NEG);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, n, h), *t_horner(y, n, h), *t_horner(z, n, h), h, step, _);
    }
    return 0;
}
