/*
 * Logistic Equation
 *
 * Example:  ./tsm-logistic-dbg 16 10 0.1 10000 .6 0 0 .1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, x0, a, h, _, *wa, *wb, *w1, *x;

    // initialize from command arguments
    assert(argc == 9);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &x0, &_, &_, &a);

    // initialize the derivative and temporary jets
    x = t_jet(n + 1);
    wa = t_jet(n);
    wb = t_jet(n);
    mpfr_set_str(_, "1.0", BASE, RND);
    w1 = t_jet_c(n, _);

    // main loop
    t_xyz_output(x0, x0, x0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(x[0], x0, RND);
        for (int k = 0; k < n; k++) {
            //  x' = Ax(1 - x)
            mpfr_mul(wa[k], a, x[k], RND);
            mpfr_sub(wb[k], w1[k], x[k], RND);
            mpfr_div_ui(x[k + 1], *t_prod(&_, wa, wb, k), k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x0, x, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x0, x0, x0, t);
    }
    return 0;
}
