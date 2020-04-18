/*
 * Logistic Equation
 *
 * Example:  ./tsm-logistic-dbg 32 10 0.1 10000 .6 0 0 .1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, a, h, _;

    // initialize from command arguments
    assert(argc == 9);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &_, &_, &a);

    // initialize the derivative and temporary jets
    mpfr_t *x = t_jet_c(n + 1, x0);

    // main loop
    t_xyz_output(x[0], x[0], x[0], t);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = Ax(1 - x)
            mpfr_fmms(_, a, x[k], a, *t_sqr(x, k), RND);
            t_next(x, _, k, POS);
        }

        // sum the series using Horner's method and advance one step
        t_horner(x, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x[0], x[0], x[0], t);
    }
    return 0;
}
