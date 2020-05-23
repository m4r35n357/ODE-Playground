/*
 * Logistic Equation
 *
 * Example:  ./tsm-logistic-dbg 9 32 10 0.1 10000 .6 0 0 .1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t x0, a, h, _;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &_, &_, &a);

    // initialize the derivative and temporary jets
    series x = t_jet_c(n + 1, x0);

    t_output(x.a[0], x.a[0], x.a[0], h, 0, _);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = Ax(1 - x)
            mpfr_fmms(_, a, x.a[k], a, *t_sqr(x, k), RND);
            t_next(x, _, k, POS);
        }

        // sum the series using Horner's method and advance one step
        t_horner(x, h);
        t_output(x.a[0], x.a[0], x.a[0], h, step, _);
    }
    return 0;
}
