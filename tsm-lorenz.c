/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-dbg 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, s, r, b, h, tmp, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 12);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    mpfr_init_set_str(x, argv[5], BASE, RND);
    mpfr_init_set_str(y, argv[6], BASE, RND);
    mpfr_init_set_str(z, argv[7], BASE, RND);
    mpfr_init_set_str(s, argv[8], BASE, RND);
    mpfr_init_set_str(r, argv[9], BASE, RND);
    mpfr_init_set_str(b, argv[10], BASE, RND);
    mpfr_init_set_str(tmp, argv[11], BASE, RND);
    mpfr_div(b, b, tmp, RND);

    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);

    // main loop
    for (long step = 1; step < nsteps + 1; step++) {
        // print a line of output
        t_line_output(t, 3, x, y, z);

        // compute the taylor coefficients
        mpfr_set(cx[0], x, RND);
        mpfr_set(cy[0], y, RND);
        mpfr_set(cz[0], z, RND);
        for (int k = 0; k < order; k++) {
            //  x' = S(y - x)
            mpfr_fmms(tmp, s, cy[k], s, cx[k], RND);
            mpfr_div_ui(cx[k + 1], tmp, k + 1, RND);
            //  y' = x(R - z) - y
            t_product(&tmp, cx, cz, k);
            mpfr_fms(tmp, r, cx[k], tmp, RND);
            mpfr_sub(tmp, tmp, cy[k], RND);
            mpfr_div_ui(cy[k + 1], tmp, k + 1, RND);
            //  z' = xy - Bz
            t_product(&tmp, cx, cy, k);
            mpfr_fms(tmp, b, cz[k], tmp, RND);
            mpfr_div_si(cz[k + 1], tmp, - (k + 1), RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x, cx, order, h);
        t_horner(&y, cy, order, h);
        t_horner(&z, cz, order, h);
        mpfr_mul_ui(t, h, step, RND);
    }
    return 0;
}
