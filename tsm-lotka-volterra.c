/*
 * Lotka-Volterra (Predator-Prey) System
 *
 * Example: ./tsm-lotka-volterra-dbg 9 32 10 .01 2000 10 10 0 1 .5 .05 .02
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t x0, y0, a, b, c, d, h, _, xy;

    // initialize from command arguments
    assert(argc == 13);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &_, &a, &b, &c, &d);
    mpfr_init(xy);

    // initialize the derivative and temporary jets
    series x = t_jet_c(n + 1, x0), y = t_jet_c(n + 1, y0);

    t_output(x.a[0], y.a[0], x.a[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            mpfr_set(xy, *t_prod(x, y, k), RND);
            //  x' = Ax - Cxy
            mpfr_fmms(_, a, x.a[k], c, xy, RND);
            t_next(x, _, k, POS);
            //  y' = Dxy - By
            mpfr_fmms(_, d, xy, b, y.a[k], RND);
            t_next(y, _, k, POS);
        }

        // sum the series using Horner's method and advance one step
        t_horner(x, h);
        t_horner(y, h);
        t_output(x.a[0], y.a[0], x.a[0], h, step);
    }
    return 0;
}
