/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-dbg 16 10 0.01 15000 1 0 0 6.7 2
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, z0, alpha, kappa, h, _;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &alpha, &kappa);
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
            //  x' = ay - kx - yz
            mpfr_fmms(_, alpha, y[k], kappa, x[k], RND);
            mpfr_sub(_, _, *t_prod(y, z, k), RND);
            mpfr_div_ui(x[k + 1], _, k + 1, RND);
            //  y' = x
            mpfr_div_ui(y[k + 1], x[k], k + 1, RND);
            //  z' = y^2 - z
            mpfr_sub(_, *t_sqr(y, k), z[k], RND);
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
