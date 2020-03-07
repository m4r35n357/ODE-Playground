/*
 * Lotka-Volterra (Predator-Prey) System
 *
 * Example: ./tsm-lotka-volterra-dbg 16 10 .01 2000 10 10 0 1 .5 .05 .02
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, a, b, c, d, h, _, xy;

    // initialize from command arguments
    assert(argc == 12);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &_, &a, &b, &c, &d);
    mpfr_init(xy);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);
    mpfr_t *y = t_jet_c(n + 1, y0);

    // main loop
    t_xyz_output(x[0], y[0], x[0], t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            t_prod(&xy, x, y, k);
            //  x' = Ax - Cxy
            mpfr_fmms(_, a, x[k], c, xy, RND);
            mpfr_div_ui(x[k + 1], _, k + 1, RND);
            //  y' = Dxy - By
            mpfr_fmms(_, d, xy, b, y[k], RND);
            mpfr_div_ui(y[k + 1], _, k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&_, x, n, h);
        t_horner(&_, y, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x[0], y[0], x[0], t);
    }
    return 0;
}
