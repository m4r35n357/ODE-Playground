/*
 * Wimol-Banlue System
 *
 * Example: ./tsm-wimol-banlue-dbg 16 8 0.1 10000 1 0 0 2.0
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, a, h, _, __;

    // initialize from command arguments
    assert(argc == 9);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a);
    mpfr_inits(_, __, NULL);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);
    mpfr_t *tx = t_jet(n);
    mpfr_t *s2x = t_jet(n);
    mpfr_t *wa = t_jet_c(n, a);

    // main loop
    t_xyz_output(x[0], y[0], z[0], t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = y - x
            mpfr_sub(_, y[k], x[k], RND);
            mpfr_div_ui(x[k + 1], _, k + 1, RND);
            //  y' = - z * tan(x)
            t_tan_sec2(tx, s2x, x, k, &_, HYP);
            mpfr_div_si(y[k + 1], *t_prod(&_, z, tx, k), - (k + 1), RND);
            //  z' = - A + xy + |y|
            mpfr_add(_, *t_prod(&_, x, y, k), *t_abs(&__, y, k), RND);
            mpfr_sub(_, _, wa[k], RND);
            mpfr_div_ui(z[k + 1], _, k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&_, x, n, h);
        t_horner(&_, y, n, h);
        t_horner(&_, z, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x[0], y[0], z[0], t);
    }
    return 0;
}
