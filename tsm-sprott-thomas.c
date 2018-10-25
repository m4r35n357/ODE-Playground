/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-dbg 16 10 0.01 30001 1 0 0 4.75 1
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, a, b, h, tmp, *wsx, *wcx, *wsy, *wcy, *wsz, *wcz, *wax, *way, *waz,
        *wsax, *wcax, *wsay, *wcay, *wsaz, *wcaz, *cx, *cy, *cz;

int main (int argc, char **argv) {
    // initialize from command arguments
    t_stepper(argc, argv, &order, &t, &h, &nsteps);
    mpfr_inits(tmp, NULL);
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
    for (long step = 1; step < nsteps + 1; step++) {
        // print a line of output
        t_line_output(t, 3, x, y, z);

        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            t_tan_sec2(wsx, wcx, cx, k, &tmp, TRIG);
            t_tan_sec2(wsy, wcy, cy, k, &tmp, TRIG);
            t_tan_sec2(wsz, wcz, cz, k, &tmp, TRIG);
            mpfr_mul(wax[k], cx[k], a, RND);
            mpfr_mul(way[k], cy[k], a, RND);
            mpfr_mul(waz[k], cz[k], a, RND);
            t_sin_cos(wsax, wcax, wax, k, &tmp, TRIG);
            t_sin_cos(wsay, wcay, way, k, &tmp, TRIG);
            t_sin_cos(wsaz, wcaz, waz, k, &tmp, TRIG);
            //  x' = sin(Ay) - Bsin(x)
            mpfr_fms(tmp, wsx[k], b, wsay[k], RND);
            mpfr_div_si(cx[k + 1], tmp, - (k + 1), RND);
            //  y' = sin(Az) - Bsin(y)
            mpfr_fms(tmp, wsy[k], b, wsaz[k], RND);
            mpfr_div_si(cy[k + 1], tmp, - (k + 1), RND);
            //  z' = sin(Ax) - Bsin(z)
            mpfr_fms(tmp, wsz[k], b, wsax[k], RND);
            mpfr_div_si(cz[k + 1], tmp, - (k + 1), RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x, cx, order, h);
        t_horner(&y, cy, order, h);
        t_horner(&z, cz, order, h);
        mpfr_mul_si(t, h, step, RND);
    }
    return 0;
}
