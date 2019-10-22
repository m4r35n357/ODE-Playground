/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-dbg 16 10 .01 50000 .1 .1 .1 .2876 .1
 *          ./tsm-rf-dbg 16 10 .01 50000 -1 0 .75 .25 .2
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, alpha, gamma, h, D3, _, wx2_1, *wa, *wb, *wc, *w1, *walpha, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 10);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x);
    t_arg(argv, 6, &y);
    t_arg(argv, 7, &z);
    t_arg(argv, 8, &alpha);
    t_arg(argv, 9, &gamma);
    mpfr_inits(_, wx2_1, NULL);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    wa = t_jet(order);
    wb = t_jet(order);
    wc = t_jet(order);
    mpfr_set_ui(_, 1, RND);
    w1 = t_jet_c(order, _);
    walpha = t_jet_c(order, alpha);
    mpfr_init_set_ui(D3, 3, RND);

    // main loop
    t_xyz_output(x, y, z, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            mpfr_sub(wx2_1, *t_sqr(&_, cx, k), w1[k], RND);
            //  x' = y(z - 1 + x^2) + Gx
            mpfr_add(wa[k], cz[k], wx2_1, RND);
            mpfr_fma(_, gamma, cx[k], *t_prod(&_, cy, wa, k), RND);
            mpfr_div_ui(cx[k + 1], _, k + 1, RND);
            //  y' = x(3z + 1 - x^2) + Gy
            mpfr_fms(wb[k], D3, cz[k], wx2_1, RND);
            mpfr_fma(_, gamma, cy[k], *t_prod(&_, cx, wb, k), RND);
            mpfr_div_ui(cy[k + 1], _, k + 1, RND);
            //  z' = -2z(A + xy)
            mpfr_add(wc[k], *t_prod(&_, cx, cy, k), walpha[k], RND);
            mpfr_mul_2ui(_, *t_prod(&_, cz, wc, k), 1, RND);
            mpfr_div_si(cz[k + 1], _, - (k + 1), RND);
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
