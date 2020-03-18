/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-dbg 16 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, sigma, rho, beta, h, _;

    // initialize from command arguments
    assert(argc == 12);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &sigma, &rho, &beta, &_);
    mpfr_div(beta, beta, _, RND);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);
    mpfr_t *z = t_jet_c(n + 1, z0);

    // main loop
    t_xyz_output(x[0], y[0], z[0], t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = S(y - x)
            mpfr_fmms(_, sigma, y[k], sigma, x[k], RND);
            mpfr_div_ui(x[k + 1], _, k + 1, RND);
            //  y' = x(R - z) - y
            mpfr_fms(_, x[k], rho, *t_prod(x, z, k), RND);
            mpfr_sub(_, _, y[k], RND);
            mpfr_div_ui(y[k + 1], _, k + 1, RND);
            //  z' = xy - Bz
            mpfr_fms(_, beta, z[k], *t_prod(x, y, k), RND);
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
