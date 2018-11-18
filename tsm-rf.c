/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-dbg 16 10 .01 100001 .1 .1 .1 .2876 .1
 *          ./tsm-rf-dbg 16 10 .01 100001 -1.0 0 .75 .25 .2
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, a, g, h, D3, _, wx2_1, *wa, *wb, *wc, *w1, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 10);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    mpfr_inits(_, wx2_1, NULL);
    mpfr_init_set_str(x, argv[5], BASE, RND);
    mpfr_init_set_str(y, argv[6], BASE, RND);
    mpfr_init_set_str(z, argv[7], BASE, RND);
    mpfr_init_set_str(a, argv[8], BASE, RND);
    mpfr_init_set_str(g, argv[9], BASE, RND);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    wa = t_jet(order);
    wb = t_jet(order);
    wc = t_jet(order);
    mpfr_set_ui(_, 1, RND);
    w1 = t_jet_constant(order, _);
    mpfr_init_set_ui(D3, 3, RND);

    // main loop
    for (long step = 1; step < nsteps + 1; step++) {
        // print a line of output
        t_line_output(t, 3, x, y, z);

        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            t_square(&_, cx, k);
            mpfr_sub(wx2_1, _, w1[k], RND);
            //  x' = y(z - 1 + x^2) + Gx
            mpfr_add(wa[k], cz[k], wx2_1, RND);
            t_product(&_, cy, wa, k);
            mpfr_fma(_, g, cx[k], _, RND);
            mpfr_div_ui(cx[k + 1], _, k + 1, RND);
            //  y' = x(3z + 1 - x^2) + Gy
            mpfr_fms(wb[k], D3, cz[k], wx2_1, RND);
            t_product(&_, cx, wb, k);
            mpfr_fma(_, g, cy[k], _, RND);
            mpfr_div_ui(cy[k + 1], _, k + 1, RND);
            //  z' = -2z(A + xy)
            t_product(&_, cx, cy, k);
            mpfr_add(wc[k], _, a, RND);
            t_product(&_, cz, wc, k);
            mpfr_mul_2ui(_, _, 1, RND);
            mpfr_div_si(cz[k + 1], _, - (k + 1), RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x, cx, order, h);
        t_horner(&y, cy, order, h);
        t_horner(&z, cz, order, h);
        mpfr_mul_ui(t, h, step, RND);
    }
    return 0;
}
