/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-dbg 16 10 .01 100001 1 0 0 1.4
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, a, h, tmp, w4x, w4y, w4z, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 9);
    // initialize from command arguments
    t_stepper(argc, argv, &order, &t, &h, &nsteps);
    mpfr_inits(tmp, w4x, w4y, w4z, NULL);
    mpfr_init_set_str(x, argv[5], BASE, RND);
    mpfr_init_set_str(y, argv[6], BASE, RND);
    mpfr_init_set_str(z, argv[7], BASE, RND);
    mpfr_init_set_str(a, argv[8], BASE, RND);

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
            mpfr_mul_2ui(w4x, cx[k], 2, RND);
            mpfr_mul_2ui(w4y, cy[k], 2, RND);
            mpfr_mul_2ui(w4z, cz[k], 2, RND);
            //  x' = - Ax - 4y - 4z - y^2
            t_square(&tmp, cy, k);
            mpfr_add(tmp, w4y, tmp, RND);
            mpfr_add(tmp, w4z, tmp, RND);
            mpfr_fma(tmp, a, cx[k], tmp, RND);
            mpfr_div_si(cx[k + 1], tmp, - (k + 1), RND);
            //  y' = - Ay - 4z - 4x - z^2
            t_square(&tmp, cz, k);
            mpfr_add(tmp, w4z, tmp, RND);
            mpfr_add(tmp, w4x, tmp, RND);
            mpfr_fma(tmp, a, cy[k], tmp, RND);
            mpfr_div_si(cy[k + 1], tmp, - (k + 1), RND);
            //  z' = - Az - 4x - 4y - x^2
            t_square(&tmp, cx, k);
            mpfr_add(tmp, w4x, tmp, RND);
            mpfr_add(tmp, w4y, tmp, RND);
            mpfr_fma(tmp, a, cz[k], tmp, RND);
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
