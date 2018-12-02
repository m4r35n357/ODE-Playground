/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-dbg 80 40 0.02 50001 1 1 0 3 2.2 1 .01
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
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
    mpfr_init_set_str(x, argv[5], BASE, RND);
    mpfr_init_set_str(y, argv[6], BASE, RND);
    mpfr_init_set_str(z, argv[7], BASE, RND);
    mpfr_init_set_str(a, argv[8], BASE, RND);
    mpfr_init_set_str(b, argv[9], BASE, RND);
    mpfr_init_set_str(g, argv[10], BASE, RND);
    mpfr_init_set_str(m, argv[11], BASE, RND);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    wa = t_jet(order);
    wb = t_jet(order);
    mpfr_init_set_ui(_, 1, RND);
    w1 = t_jet_c(order, _);

    // main loop
    for (long step = 1; step < nsteps + 1; step++) {
        // print a line of output
        t_line_output(t, 3, x, y, z);

        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            //  x' = Ax(1 - y) - Bz
            mpfr_sub(wa[k], w1[k], cy[k], RND);
            t_product(&_, cx, wa, k);
            mpfr_fmms(_, a, _, b, cz[k], RND);
            mpfr_div_ui(cx[k + 1], _, k + 1, RND);
            //  y' = - Gy(1 - x^2)
            t_square(&_, cx, k);
            mpfr_sub(wb[k], w1[k], _, RND);
            t_product(&_, cy, wb, k);
            mpfr_mul(_, _, g, RND);
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
    }
    return 0;
}
