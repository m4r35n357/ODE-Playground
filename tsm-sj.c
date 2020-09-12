/*
 * Sprott-Jafari System
 *
 * Example: ./tsm-sj-dbg NA NA 10 .01 10000 0 3.9 .7 8.888 4
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    real a, h;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_jet(n + 1), y = t_jet(n + 1), z = t_jet(n + 1), w_b = t_jet(n);
    t_args(argv, argc, x, y, z, &a, w_b);

    // main loop
    t_xyz_output(x[0], y[0], z[0], 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            x[k + 1] = y[k] / (k + 1);
            y[k + 1] = - x[k] + t_prod(y, z, k) / (k + 1);
            z[k + 1] = (z[k] + a * t_sqr(x, k) - t_sqr(y, k) - w_b[k]) / (k + 1);
        }

        // sum the series using Horner's method and advance one step
        t_xyz_output(t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step);
    }
    return 0;
}
