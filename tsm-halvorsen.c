/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-dbg NA NA 10 .01 10000 1 0 0 1.4
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
    assert(argc == 10);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_jet(n + 1), y = t_jet(n + 1), z = t_jet(n + 1);
    t_args(argv, argc, x, y, z, &a);

    // main loop
    t_xyz_output(x[0], y[0], z[0], 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            x[k + 1] = - (a * x[k] + 4.0 * (y[k] + z[k]) + t_sqr(y, k)) / (k + 1);
            y[k + 1] = - (a * y[k] + 4.0 * (z[k] + x[k]) + t_sqr(z, k)) / (k + 1);
            z[k + 1] = - (a * z[k] + 4.0 * (x[k] + y[k]) + t_sqr(x, k)) / (k + 1);
        }

        // sum the series using Horner's method and advance one step
        t_xyz_output(t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step);
    }
    return 0;
}
