/*
 * Genesio-Tesi System
 *
 * Example: ./tsm-genesio-tesi-dbg NA NA 10 0.01 50000 .1 .1 .1 .44 1.1 1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    long double a, b, c, h;

    // initialize from command arguments
    assert(argc == 12);
    t_stepper(argv, &n, &h, &nsteps);
    long double *x = t_jet(n + 1), *y = t_jet(n + 1), *z = t_jet(n + 1);
    t_args(argv, argc, x, y, z, &a, &b, &c);

    // main loop
    t_xyz_output(x[0], y[0], z[0], 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            x[k + 1] = y[k] / (k + 1);
            y[k + 1] = z[k] / (k + 1);
            z[k + 1] = (t_sqr(x, k) - c * x[k] - b * y[k] - a * z[k]) / (k + 1);
        }

        // sum the series using Horner's method and advance one step
        t_xyz_output(t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step);
    }
    return 0;
}
