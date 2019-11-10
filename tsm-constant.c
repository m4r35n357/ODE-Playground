/*
 * Proportional
 *
 * Example: ./tsm-constant-dbg 16 10 0.1 10001 10 0 0 -.05
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, a, h, _, *cx;

int main (int argc, char **argv) {
    assert(argc == 9);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x);
    t_arg(argv, 8, &a);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);

    // main loop
    t_xyz_output(x, x, x, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        for (int k = 0; k < order; k++) {
            //  x' = Ax
            mpfr_mul(_, cx[k], a, RND);
            mpfr_div_ui(cx[k + 1], _, k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x, cx, order, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x, x, x, t);
    }
    return 0;
}