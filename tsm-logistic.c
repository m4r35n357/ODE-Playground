/*
 * Logistic Equation
 *
 * Example:  ./tsm-logistic-dbg 9 32 10 0.005 10000 .001 0 0 .5
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
    t_stepper(argv, &n, &h, &nsteps);
    mpfr_inits(a, _, NULL);
    series x = t_series(n + 1);
    t_args(argv, argc, x.jet, &_, &_, &a);

    t_output(x.jet[0], _, _, h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = Ax(1 - x)
            mpfr_fmms(_, a, x.jet[k], a, *t_sqr(x, k), RND);
            t_next(x, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), _, _, h, step);
    }
    return 0;
}
