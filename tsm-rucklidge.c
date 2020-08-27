/*
 * Rucklidge Attractor
 *
 * Example: ./tsm-rucklidge-dbg NA NA 10 0.01 15000 1 0 0 6.7 2
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    long double alpha, kappa, h;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &h, &nsteps);
    long double *x = t_jet(n + 1), *y = t_jet(n + 1), *z = t_jet(n + 1);
    t_args(argv, argc, x, y, z, &alpha, &kappa);

    // main loop
    t_xyz_output(x[0], y[0], z[0], 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            x[k + 1] = (alpha * y[k] - kappa * x[k] - t_prod(y, z, k)) / (k + 1);
            y[k + 1] = x[k] / (k + 1);
            z[k + 1] = (t_sqr(y, k) - z[k]) / (k + 1);
        }

        // sum the series using Horner's method and advance one step
        t_xyz_output(t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step);
    }
    return 0;
}
