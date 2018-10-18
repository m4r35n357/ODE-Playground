/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-dbg 16 10 0.01 10001 1 0 0 6.7 2
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, alpha, kappa, h, tmp, tmp2, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 10);
    // initialize from command arguments
    t_stepper(argc, argv, &order, &t, &h, &nsteps);
    mpfr_inits(tmp, tmp2, NULL);
    mpfr_init_set_str(x, argv[5], BASE, RND);
    mpfr_init_set_str(y, argv[6], BASE, RND);
    mpfr_init_set_str(z, argv[7], BASE, RND);
    mpfr_init_set_str(alpha, argv[8], BASE, RND);
    mpfr_init_set_str(kappa, argv[9], BASE, RND);

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
            //  x' = ay - kx - yz
            t_product(&tmp, cy, cz, k);
            mpfr_fmms(tmp2, alpha, cy[k], kappa, cx[k], RND);
            mpfr_sub(tmp, tmp2, tmp, RND);
            mpfr_div_si(cx[k + 1], tmp, k + 1, RND);
            //  y' = x
            mpfr_div_si(cy[k + 1], cx[k], k + 1, RND);
            //  z' = y^2 - z
            t_square(&tmp, cy, k);
            mpfr_sub(tmp, tmp, cz[k], RND);
            mpfr_div_si(cz[k + 1], tmp, k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&x, cx, order, h);
        t_horner(&y, cy, order, h);
        t_horner(&z, cz, order, h);
        mpfr_mul_si(t, h, step, RND);
    }
    return 0;
}
