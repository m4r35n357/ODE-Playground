/*
 * Wimol-Banlue System
 *
 * Example: ./tsm-wimol-banlue-dbg 16 8 0.1 10001 1 0 0 2.0
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, a, b, c, h, _, __, *wtx, *w__, *w_a, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 9);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    mpfr_inits(_, __, NULL);
    mpfr_init_set_str(x, argv[5], BASE, RND);
    mpfr_init_set_str(y, argv[6], BASE, RND);
    mpfr_init_set_str(z, argv[7], BASE, RND);
    mpfr_init_set_str(a, argv[8], BASE, RND);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    wtx = t_jet(order);
    w__ = t_jet(order);
    w_a = t_jet_c(order, a);

    // main loop
    t_xyz_output(x, y, z, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            t_tan_sec2(wtx, w__, cx, k, &_, HYP);
            //  x' = y - x
            mpfr_sub(_, cy[k], cx[k], RND);
            mpfr_div_ui(cx[k + 1], _, k + 1, RND);
            //  y' = - z * tan(x)
            t_prod(&_, cz, wtx, k);
            mpfr_div_si(cy[k + 1], _, - (k + 1), RND);
            //  z' = - A + x * y + |y|
            t_prod(&_, cx, cy, k);
            t_abs(&__, cy, k);
            mpfr_add(_, _, __, RND);
            mpfr_sub(_, _, w_a[k], RND);
            mpfr_div_ui(cz[k + 1], _, k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x, cx, order, h);
        t_horner(&y, cy, order, h);
        t_horner(&z, cz, order, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x, y, z, t);
    }
    return 0;
}
