/*
 * Proportional
 *
 * Example: ./tsm-constant-dbg NA NA 10 0.1 10000 10 0 0 -.05
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    long double a, h, _;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &h, &nsteps);
    long double *x = t_jet(n + 1);
    t_args(argv, argc, x, &_, &_, &a);

    // main loop
    t_xyz_output(x[0], 0.0, 0.0, 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            x[k + 1] = a * x[k] / (k + 1);
        }

        // sum the series using Horner's method and advance one step
        t_xyz_output(t_horner(x, n, h), 0.0, 0.0, h * step);
    }
    return 0;
}
