/*
 * Yu-Wang System
 *
 * Example: ./tsm-yu-wang-dbg 9 32 10 .001 50000 1 0 0 10 40 2 2.5
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t x0, y0, z0, a, b, c, d, h, _;

    // initialize from command arguments
    assert(argc == 13);
    t_stepper(argv, &n, &h, &nsteps);
    t_args(argv, argc, &x0, &y0, &z0, &a, &b, &c, &d);
    mpfr_init(_);

    // initialize the derivative and temporary jets
    series x = t_jet_c(n + 1, x0), y = t_jet_c(n + 1, y0), z = t_jet_c(n + 1, z0);
    series xy = t_jet(n), e_xy = t_jet(n);

    t_output(x.a[0], y.a[0], z.a[0], h, 0, _);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = A(y - x)
            mpfr_fmms(_, a, y.a[k], a, x.a[k], RND);
            t_next(x, _, k, POS);
            //  y' = Bx - cxz
            mpfr_fmms(_, b, x.a[k], c, *t_prod(x, z, k), RND);
            t_next(y, _, k, POS);
            //  z' = e^(xy) - Dz
            mpfr_set(xy.a[k], *t_prod(x, y, k), RND);
            mpfr_fms(_, d, z.a[k], *t_exp(e_xy, xy, k), RND);
            t_next(z, _, k, NEG);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step, _);
    }
    return 0;
}
