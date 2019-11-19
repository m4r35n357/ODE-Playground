/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-dbg 16 10 0.01 15000 1 0 0 6.7 2
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long n, nsteps;
mpfr_t t, x0, y0, z0, alpha, kappa, h, _, *x, *y, *z;

int main (int argc, char **argv) {
    assert(argc == 10);
    // initialize from command arguments
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 6, &y0);
    t_arg(argv, 7, &z0);
    t_arg(argv, 8, &alpha);
    t_arg(argv, 9, &kappa);
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
            //  x' = ay - kx - yz
            mpfr_fms(_, alpha, y[k], *t_prod(&_, y, z, k), RND);
            mpfr_fms(_, kappa, x[k], _, RND);
            mpfr_div_si(x[k + 1], _, - (k + 1), RND);
            //  y' = x
            mpfr_div_ui(y[k + 1], x[k], k + 1, RND);
            //  z' = y^2 - z
            mpfr_sub(_, *t_sqr(&_, y, k), z[k], RND);
            mpfr_div_ui(z[k + 1], _, k + 1, RND);
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
