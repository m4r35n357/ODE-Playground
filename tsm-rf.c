/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-dbg NA NA 10 .01 50000 .05 -.05 .3 .28713 .1
 *          ./tsm-rf-dbg NA NA 16 .01 50000 .05 -.05 .3 .105 .1 | ./plotPi3d.py
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    real gamma, h, x2_1;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &h, &nsteps);
    series x = t_jet(n + 1), y = t_jet(n + 1), z = t_jet(n + 1), alpha = t_jet(n);
    t_args(argv, argc, x, y, z, alpha, &gamma);
    series a = t_jet(n), b = t_jet(n), c = t_jet(n), w1 = t_jet_c(n, 1.0);

    // main loop
    t_xyz_output(x[0], y[0], z[0], 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            x2_1 = t_sqr(x, k) - w1[k];
            a[k] = z[k] + x2_1;
            b[k] = 3.0 * z[k] - x2_1;
            c[k] = alpha[k] + t_prod(x, y, k);
            x[k + 1] = (t_prod(y, a, k) + gamma * x[k]) / (k + 1);
            y[k + 1] = (t_prod(x, b, k) + gamma * y[k]) / (k + 1);
            z[k + 1] = - 2.0 * t_prod(z, c, k) / (k + 1);
        }

        // sum the series using Horner's method and advance one step
        t_xyz_output(t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step);
    }
    return 0;
}
