/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-dbg NA NA 16 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    real sigma, rho, beta, h, _;

    // initialize from command arguments
    assert(argc == 13);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_jet(n + 1), y = t_jet(n + 1), z = t_jet(n + 1);
    t_args(argv, argc, x, y, z, &sigma, &rho, &beta, &_);
    beta /= _;

    // main loop
    t_xyz_output(x[0], y[0], z[0], 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            x[k + 1] = sigma * (y[k] - x[k]) / (k + 1);
            y[k + 1] = (rho * x[k] - y[k] - t_prod(x, z, k)) / (k + 1);
            z[k + 1] = (t_prod(x, y, k) - beta * z[k]) / (k + 1);
        }

        // sum the series using Horner's method and advance one step
        t_xyz_output(t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step);
    }
    return 0;
}
