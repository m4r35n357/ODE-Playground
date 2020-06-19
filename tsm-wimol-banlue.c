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
    mpfr_init(_);
    series x = t_series(n + 1), y = t_series(n + 1), z = t_series(n + 1);
    series wa = t_series(n);
    t_args(argv, argc,x.jet, y.jet, z.jet, wa.jet);
    series tx = t_series(n), s2x = t_series(n);

    t_output(x.jet[0], y.jet[0], z.jet[0], h, 0);
    for (long step = 1; step <= nsteps; step++) {
        // build the jet of taylor coefficients
        for (int k = 0; k < n; k++) {
            //  x' = y - x
            mpfr_sub(_, y.jet[k], x.jet[k], RND);
            t_next(x, _, k, POS);
            //  y' = - z * tan(x)
            t_tan_sec2(tx, s2x, x, k, HYP);
            t_next(y, *t_prod(z, tx, k), k, NEG);
            //  z' = - A + xy + |y|
            mpfr_add(_, *t_prod(x, y, k), *t_abs(y, k), RND);
            mpfr_sub(_, _, wa.jet[k], RND);
            t_next(z, _, k, POS);
        }
        // sum the series using Horner's method and advance one step
        t_output(*t_horner(x, h), *t_horner(y, h), *t_horner(z, h), h, step);
    }
    return 0;
}
