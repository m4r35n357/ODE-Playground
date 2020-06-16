/*
 * Exponential growth & decay
 *
 * Example: ./tsm-constant-dbg 9 32 10 0.1 10000 10 0 0 -.05
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t a, h, _;

    // initialize from command arguments
    assert(argc == 10);
    mpfr_inits(a, _, NULL);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_series(n + 1);
    t_args(argv, argc, x.jet, &_, &_, &a);

    t_output(x.jet[0], x.jet[0], x.jet[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = Ax
            mpfr_mul(_, x.jet[k], a, RND);
            t_next(x, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_horner(x, h);
        t_output(x.jet[0], x.jet[0], x.jet[0], h, step);
    }
    return 0;
}
