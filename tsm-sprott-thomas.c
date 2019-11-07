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

long n, nsteps;
mpfr_t t, x0, y0, z0, a, b, h, _, *wsx, *wcx, *wsy, *wcy, *wsz, *wcz, *wax, *way, *waz,
        *wsax, *wcax, *wsay, *wcay, *wsaz, *wcaz, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 10);
    // initialize from command arguments
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_arg(argv, 5, &x0);
    t_arg(argv, 6, &y0);
    t_arg(argv, 7, &z0);
    t_arg(argv, 8, &a);
    t_arg(argv, 9, &b);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    cx = t_jet(n + 1);
    cy = t_jet(n + 1);
    cz = t_jet(n + 1);
    wax = t_jet(n);
    way = t_jet(n);
    waz = t_jet(n);
    wsax = t_jet(n);
    wsay = t_jet(n);
    wsaz = t_jet(n);
    wcax = t_jet(n);
    wcay = t_jet(n);
    wcaz = t_jet(n);
    wsx = t_jet(n);
    wsy = t_jet(n);
    wsz = t_jet(n);
    wcx = t_jet(n);
    wcy = t_jet(n);
    wcz = t_jet(n);

    // main loop
    t_xyz_output(x0, y0, z0, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x0, RND);
        mpfr_set(cy[0], y0, RND);
        mpfr_set(cz[0], z0, RND);
        for (int k = 0; k < n; k++) {
            mpfr_mul(wax[k], cx[k], a, RND);
            mpfr_mul(way[k], cy[k], a, RND);
            mpfr_mul(waz[k], cz[k], a, RND);
            //  x' = sin(Ay) - Bsin(x)
            mpfr_fms(_, t_tan_sec2(wsx, wcx, cx, k, &_, TRIG).a[k], b, t_sin_cos(wsay, wcay, way, k, &_, TRIG).a[k], RND);
            mpfr_div_si(cx[k + 1], _, - (k + 1), RND);
            //  y' = sin(Az) - Bsin(y)
            mpfr_fms(_, t_tan_sec2(wsy, wcy, cy, k, &_, TRIG).a[k], b, t_sin_cos(wsaz, wcaz, waz, k, &_, TRIG).a[k], RND);
            mpfr_div_si(cy[k + 1], _, - (k + 1), RND);
            //  z' = sin(Ax) - Bsin(z)
            mpfr_fms(_, t_tan_sec2(wsz, wcz, cz, k, &_, TRIG).a[k], b, t_sin_cos(wsax, wcax, wax, k, &_, TRIG).a[k], RND);
            mpfr_div_si(cz[k + 1], _, - (k + 1), RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x0, cx, n, h);
        t_horner(&y0, cy, n, h);
        t_horner(&z0, cz, n, h);
        mpfr_mul_ui(t, h, step, RND);
        t_xyz_output(x0, y0, z0, t);
    }
    return 0;
}
