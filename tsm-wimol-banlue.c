/*
 * Wimol-Banlue System
 *
 * Example: ./tsm-wimol-banlue-dbg NA NA 8 0.1 10000 1 0 0 2.0
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    real h;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_jet(n + 1), y = t_jet(n + 1), z = t_jet(n + 1), wa = t_jet(n);
    t_args(argv, argc, x, y, z, wa);
    series tx = t_jet(n), s2x = t_jet(n);

    // main loop
    t_xyz_output(x[0], y[0], z[0], 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            tx[k] = t_tan_sec2(tx, s2x, x, k, HYP).a;
            x[k + 1] = (y[k] - x[k]) / (k + 1);
            y[k + 1] = - t_prod(z, tx, k) / (k + 1);
            z[k + 1] = (- wa[k] + t_prod(x, y, k) + t_abs(y, k)) / (k + 1);
        }

        // sum the series using Horner's method and advance one step
        t_xyz_output(t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step);
    }
    return 0;
}
