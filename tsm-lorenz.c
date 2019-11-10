/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
 *
 *          ./tsm-lorenz-static 16 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 >/tmp/dataA
 *          ./tsm-lorenz-static 16 10 .005 20000 -15.8 -17.48 35.64 10 28 8 3 | sed -n '1~2p' >/tmp/dataB
 *          ./compare.py /tmp/dataA /tmp/dataB 3
 *          ./divergence.py /tmp/dataA /tmp/dataB 3 1.0
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long n, nsteps;
mpfr_t t, x0, y0, z0, sigma, rho, beta, h, _, *x, *y, *z;

int main (int argc, char **argv) {
    assert(argc == 12);
    // initialize from command arguments
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 6, &y0);
    t_arg(argv, 7, &z0);
    t_arg(argv, 8, &sigma);
    t_arg(argv, 9, &rho);
    t_arg(argv, 10, &beta);
    t_arg(argv, 11, &_);
    mpfr_div(beta, beta, _, RND);

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
            //  x' = S(y - x)
            mpfr_fmms(_, sigma, y[k], sigma, x[k], RND);
            mpfr_div_ui(x[k + 1], _, k + 1, RND);
            //  y' = x(R - z) - y
            mpfr_fms(_, x[k], rho, *t_prod(&_, x, z, k), RND);
            mpfr_sub(_, _, y[k], RND);
            mpfr_div_ui(y[k + 1], _, k + 1, RND);
            //  z' = xy - Bz
            mpfr_fms(_, beta, z[k], *t_prod(&_, x, y, k), RND);
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