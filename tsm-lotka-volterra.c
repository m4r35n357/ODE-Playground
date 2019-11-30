/*
 * Lotka-Volterra (Predator-Prey) System
 *
 * Example: ./tsm-lotka-volterra-dbg 16 10 .01 2000 10 10 0 1 .5 .05 .02
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, y0, a, b, c, d, h, _, xy, *x, *y;

    // initialize from command arguments
    assert(argc == 12);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 6, &y0);
    t_arg(argv, 8, &a);
    t_arg(argv, 9, &b);
    t_arg(argv, 10, &c);
    t_arg(argv, 11, &d);
    mpfr_inits(_, xy, NULL);

    // initialize the derivative and temporary jets
    x = t_jet(n + 1);
    y = t_jet(n + 1);

    // main loop
    t_xyz_output(x0, y0, x0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        mpfr_set(y[0], y0, RND);
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
        t_horner(&x0, x, n, h);
        t_horner(&y0, y, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x0, y0, x0, t);
    }
    return 0;
}
