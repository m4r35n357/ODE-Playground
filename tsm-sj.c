/*
 * Sprott-Jafari System
 *
 * Example: ./tsm-sj-dbg 32 10 .01 10000 0 3.9 .7 8.888 4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, a, b, h, _, __;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b);
    mpfr_inits(_, __, NULL);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);
    mpfr_t *w_b = t_jet_c(n, b);

    // main loop
    t_xyz_output(x[0], y[0], z[0], t);
    for (long step = 1; step <= nsteps; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = y
            mpfr_div_ui(x[k + 1], y[k], k + 1, RND);
            //  y' = yz - x
            mpfr_sub(_, *t_prod(y, z, k), x[k], RND);
            mpfr_div_ui(y[k + 1], _, k + 1, RND);
            //  z' = z - ax^2 - y^2 - b
            mpfr_sub(__, z[k], w_b[k], RND);
            mpfr_sub(__, __, *t_sqr(y, k), RND);
            mpfr_fma(_, *t_sqr(x, k), a, __, RND);
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
