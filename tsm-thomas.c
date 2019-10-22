/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-dbg 16 10 0.1 30001 1 0 0 .19
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, b, h, _, *wsy, *wcy, *wsz, *wcz, *wsx, *wcx, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 9);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x);
    t_arg(argv, 6, &y);
    t_arg(argv, 7, &z);
    t_arg(argv, 8, &b);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    wsx = t_jet(order);
    wsy = t_jet(order);
    wsz = t_jet(order);
    wcx = t_jet(order);
    wcy = t_jet(order);
    wcz = t_jet(order);

    // main loop
    t_xyz_output(x, y, z, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            t_sin_cos(wsx, wcx, cx, k, &_, TRIG);
            t_sin_cos(wsy, wcy, cy, k, &_, TRIG);
            t_sin_cos(wsz, wcz, cz, k, &_, TRIG);
            //  x' = sin(y) - Bx
            mpfr_fms(_, cx[k], b, wsy[k], RND);
            mpfr_div_si(cx[k + 1], _, - (k + 1), RND);
            //  y' = sin(z) - By
            mpfr_fms(_, cy[k], b, wsz[k], RND);
            mpfr_div_si(cy[k + 1], _, - (k + 1), RND);
            //  z' = sin(x) - Bz
            mpfr_fms(_, cz[k], b, wsx[k], RND);
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
