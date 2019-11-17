/*
 * Logistic Equation
 *
 * Example:  ./tsm-logistic-dbg 16 10 0.1 10000 .6 0 0 .1
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, a, h, _, *wa, *wb, *w1, *cx;

int main (int argc, char **argv) {
    assert(argc == 9);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x);
    t_arg(argv, 8, &a);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    wa = t_jet(order);
    wb = t_jet(order);
    mpfr_set_str(_, "1.0", BASE, RND);
    w1 = t_jet_c(order, _);

    // main loop
    t_xyz_output(x, x, x, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        for (int k = 0; k < order; k++) {
            //  x' = Ax(1 - x)
            mpfr_mul(wa[k], a, cx[k], RND);
            mpfr_sub(wb[k], w1[k], cx[k], RND);
            mpfr_div_ui(cx[k + 1], *t_prod(&_, wa, wb, k), k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x, cx, order, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x, x, x, t);
    }
    return 0;
}
