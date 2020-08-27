/*
 * Lotka-Volterra (Predator-Prey) System
 *
 * Example: ./tsm-lotka-volterra-dbg NA NA 10 .01 2000 10 10 0 1 .5 .05 .02
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    long double a, b, c, d, h, xy, _;

    // initialize from command arguments
    assert(argc == 13);
    t_stepper(argv, &n, &h, &nsteps);
    long double *x = t_jet(n + 1), *y = t_jet(n + 1);
    t_args(argv, argc, x, y, &_, &a, &b, &c, &d);

    // main loop
    t_xyz_output(x[0], y[0], 0.0, 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            xy = t_prod(x, y, k);
            x[k + 1] = (a * x[k] - c * xy) / (k + 1);
            y[k + 1] = (d * xy - b * y[k]) / (k + 1);
        }

        // sum the series using Horner's method and advance one step
        t_xyz_output(t_horner(x, n, h), t_horner(y, n, h), 0.0, h * step);
    }
    return 0;
}
