/*
 * Proportional
 *
 * Example: ./tsm-constant-dbg 16 10 0.1 10000 10 0 0 -.05
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, a, h, _, *x;

    // initialize from command arguments
    assert(argc == 9);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 8, &a);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    x = t_jet(n + 1);

    // main loop
    t_xyz_output(x0, x0, x0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        for (int k = 0; k < n; k++) {
            //  x' = Ax
            mpfr_mul(_, x[k], a, RND);
            mpfr_div_ui(x[k + 1], _, k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x0, x, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x0, x0, x0, t);
    }
    return 0;
}
