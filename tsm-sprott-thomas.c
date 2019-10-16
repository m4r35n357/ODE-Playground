/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-dbg 16 10 0.01 30001 1 0 0 4.75 1
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, a, b, h, _, *wsx, *wcx, *wsy, *wcy, *wsz, *wcz, *wax, *way, *waz,
        *wsax, *wcax, *wsay, *wcay, *wsaz, *wcaz, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 10);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    mpfr_inits(_, NULL);
    mpfr_init_set_str(x, argv[5], BASE, RND);
    mpfr_init_set_str(y, argv[6], BASE, RND);
    mpfr_init_set_str(z, argv[7], BASE, RND);
    mpfr_init_set_str(a, argv[8], BASE, RND);
    mpfr_init_set_str(b, argv[9], BASE, RND);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    wax = t_jet(order);
    way = t_jet(order);
    waz = t_jet(order);
    wsax = t_jet(order);
    wsay = t_jet(order);
    wsaz = t_jet(order);
    wcax = t_jet(order);
    wcay = t_jet(order);
    wcaz = t_jet(order);
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
            t_tan_sec2(wsx, wcx, cx, k, &_, TRIG);
            t_tan_sec2(wsy, wcy, cy, k, &_, TRIG);
            t_tan_sec2(wsz, wcz, cz, k, &_, TRIG);
            mpfr_mul(wax[k], cx[k], a, RND);
            mpfr_mul(way[k], cy[k], a, RND);
            mpfr_mul(waz[k], cz[k], a, RND);
            t_sin_cos(wsax, wcax, wax, k, &_, TRIG);
            t_sin_cos(wsay, wcay, way, k, &_, TRIG);
            t_sin_cos(wsaz, wcaz, waz, k, &_, TRIG);
            //  x' = sin(Ay) - Bsin(x)
            mpfr_fms(_, wsx[k], b, wsay[k], RND);
            mpfr_div_si(cx[k + 1], _, - (k + 1), RND);
            //  y' = sin(Az) - Bsin(y)
            mpfr_fms(_, wsy[k], b, wsaz[k], RND);
            mpfr_div_si(cy[k + 1], _, - (k + 1), RND);
            //  z' = sin(Ax) - Bsin(z)
            mpfr_fms(_, wsz[k], b, wsax[k], RND);
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
