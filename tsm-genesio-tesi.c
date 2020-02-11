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
    mpfr_t t, x0, y0, z0, a, b, c, h, _, *x, *y, *z;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b, &c);
    mpfr_init(_);

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
            //  x' = y
            mpfr_div_ui(x[k + 1], y[k], k + 1, RND);
            //  y' = z
            mpfr_div_ui(y[k + 1], z[k], k + 1, RND);
            //  z' = x^2 - Cx - By -Az
            mpfr_fms(_, c, x[k], *t_sqr(&_, x, k), RND);
            mpfr_fma(_, b, y[k], _, RND);
            mpfr_fma(_, a, z[k], _, RND);
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
