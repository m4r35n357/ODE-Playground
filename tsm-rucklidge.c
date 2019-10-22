/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-dbg 16 10 0.01 10001 1 0 0 6.7 2
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

long order, nsteps;
mpfr_t t, x, y, z, alpha, kappa, h, _, __, *cx, *cy, *cz;

int main (int argc, char **argv) {
    assert(argc == 10);
    // initialize from command arguments
    t_stepper(argv, &order, &t, &h, &nsteps);
    t_arg(argv, 5, &x);
    t_arg(argv, 6, &y);
    t_arg(argv, 7, &z);
    t_arg(argv, 8, &alpha);
    t_arg(argv, 9, &kappa);
    mpfr_inits(_, __, NULL);

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
            //  x' = ay - kx - yz
            mpfr_fmms(__, alpha, cy[k], kappa, cx[k], RND);
            mpfr_sub(_, __, *t_prod(&_, cy, cz, k), RND);
            mpfr_div_ui(cx[k + 1], _, k + 1, RND);
            //  y' = x
            mpfr_div_ui(cy[k + 1], cx[k], k + 1, RND);
            //  z' = y^2 - z
            mpfr_sub(_, *t_sqr(&_, cy, k), cz[k], RND);
            mpfr_div_ui(cz[k + 1], _, k + 1, RND);
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
