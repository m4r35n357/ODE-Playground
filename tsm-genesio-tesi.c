/*
 * Genesio-Tesi System
 *
 * Example: ./tsm-genesio-tesi-dbg 16 10 0.01 150000 .1 .1 .1 .44 1.1 1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, a, b, c, h, _;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b, &c);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);

    // main loop
    t_xyz_output(x[0], y[0], z[0], t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = y
            mpfr_div_ui(x[k + 1], y[k], k + 1, RND);
            //  y' = z
            mpfr_div_ui(y[k + 1], z[k], k + 1, RND);
            //  z' = x^2 - Cx - By -Az
            mpfr_fms(_, c, x[k], *t_sqr(x, k), RND);
            mpfr_fma(_, b, y[k], _, RND);
            mpfr_fma(_, a, z[k], _, RND);
            mpfr_div_si(z[k + 1], _, - (k + 1), RND);
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
