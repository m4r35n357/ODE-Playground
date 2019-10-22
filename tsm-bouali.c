/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-dbg 80 40 0.02 50001 1 1 0 3 2.2 1 .01
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, a, b, g, m, h, _, *wa, *wb, *w1, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 12);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x);
    t_arg(argv, 6, &y);
    t_arg(argv, 7, &z);
    t_arg(argv, 8, &a);
    t_arg(argv, 9, &b);
    t_arg(argv, 10, &g);
    t_arg(argv, 11, &m);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    wa = t_jet(order);
    wb = t_jet(order);
    mpfr_init_set_ui(_, 1, RND);
    w1 = t_jet_c(order, _);

    // main loop
    t_xyz_output(x, y, z, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            //  x' = Ax(1 - y) - Bz
            mpfr_sub(wa[k], w1[k], cy[k], RND);
            mpfr_fmms(_, a, *t_prod(&_, cx, wa, k), b, cz[k], RND);
            mpfr_div_ui(cx[k + 1], _, k + 1, RND);
            //  y' = - Gy(1 - x^2)
            mpfr_sub(wb[k], w1[k], *t_sqr(&_, cx, k), RND);
            mpfr_mul(_, *t_prod(&_, cy, wb, k), g, RND);
            mpfr_div_si(cy[k + 1], _, - (k + 1), RND);
            //  z' = Mx
            mpfr_mul(_, m, cx[k], RND);
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
