/*
 * Wimol-Banlue System
 *
 * Example: ./tsm-wimol-banlue-dbg 9 32 8 0.1 10000 1 0 0 2.0
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t h, _;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_jet(n + 1), y = t_jet(n + 1), z = t_jet(n + 1);
    series wa = t_jet(n);
    t_args(argv, argc,x.a, y.a, z.a, wa.a);
    series tx = t_jet(n), s2x = t_jet(n);
    mpfr_init(_);

    t_output(x.a[0], y.a[0], z.a[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = y - x
            mpfr_sub(_, y.a[k], x.a[k], RND);
            t_next(x, _, k, POS);
            //  y' = - z * tan(x)
            t_tan_sec2(tx, s2x, x, k, HYP);
            t_next(y, *t_prod(z, tx, k), k, NEG);
            //  z' = - A + xy + |y|
            mpfr_add(_, *t_prod(x, y, k), *t_abs(y, k), RND);
            mpfr_sub(_, _, wa.a[k], RND);
            t_next(z, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
