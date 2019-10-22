/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
 *
 *          ./tsm-lorenz-static 16 10 .01 10000 -15.8 -17.48 35.64 10 28 8 3 >/tmp/dataA
 *          ./tsm-lorenz-static 16 10 .005 20000 -15.8 -17.48 35.64 10 28 8 3 | sed -n '1~2p' >/tmp/dataB
 *          ./compare.py /tmp/dataA /tmp/dataB 3
 *          ./divergence.py /tmp/dataA /tmp/dataB 3 1.0
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, s, r, b, h, _, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 12);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x);
    t_arg(argv, 6, &y);
    t_arg(argv, 7, &z);
    t_arg(argv, 8, &s);
    t_arg(argv, 9, &r);
    t_arg(argv, 10, &b);
    t_arg(argv, 11, &_);
    mpfr_div(b, b, _, RND);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);

    // main loop
    t_xyz_output(x, y, z, t);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            //  x' = S(y - x)
            mpfr_fmms(_, s, cy[k], s, cx[k], RND);
            mpfr_div_ui(cx[k + 1], _, k + 1, RND);
            //  y' = x(R - z) - y
            mpfr_fms(_, cx[k], r, *t_prod(&_, cx, cz, k), RND);
            mpfr_sub(_, _, cy[k], RND);
            mpfr_div_ui(cy[k + 1], _, k + 1, RND);
            //  z' = xy - Bz
            mpfr_fms(_, b, cz[k], *t_prod(&_, cx, cy, k), RND);
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
