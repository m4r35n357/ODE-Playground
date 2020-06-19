/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-dbg 9 32 10 .01 50000 .05 -.05 .3 .28713 .1
 *          ./tsm-rf-dbg 9 32 16 .01 50000 .05 -.05 .3 .105 .1 | ./plot3d.py
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t gamma, h, d3, _, x2_1;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &h, &nsteps);
    mpfr_inits(gamma, d3, x2_1, _, NULL);
    series x = t_series(n + 1), y = t_series(n + 1), z = t_series(n + 1), alpha = t_series(n);
    t_args(argv, argc, x.jet, y.jet, z.jet, alpha.jet, &gamma);
    series a = t_series(n), b = t_series(n), c = t_series(n), w1 = t_series(n);
    mpfr_set_ui(w1.jet[0], 1, RND);
    mpfr_set_ui(d3, 3, RND);

    t_output(x.jet[0], y.jet[0], z.jet[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            mpfr_sub(x2_1, *t_sqr(x, k), w1.jet[k], RND);
            //  x' = y(z - 1 + x^2) + Gx
            mpfr_add(a.jet[k], z.jet[k], x2_1, RND);
            mpfr_fma(_, gamma, x.jet[k], *t_prod(y, a, k), RND);
            t_next(x, _, k, POS);
            //  y' = x(3z + 1 - x^2) + Gy
            mpfr_fms(b.jet[k], d3, z.jet[k], x2_1, RND);
            mpfr_fma(_, gamma, y.jet[k], *t_prod(x, b, k), RND);
            t_next(y, _, k, POS);
            //  z' = -2z(A + xy)
            mpfr_add(c.jet[k], *t_prod(x, y, k), alpha.jet[k], RND);
            mpfr_mul_2ui(_, *t_prod(z, c, k), 1, RND);
            t_next(z, _, k, NEG);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
