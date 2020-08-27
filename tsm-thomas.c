/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-dbg NA NA 10 0.1 30000 1 0 0 .19
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    long double b, h;

    // initialize from command arguments
    assert(argc == 10);
    t_stepper(argv, &n, &h, &nsteps);
    long double *x = t_jet(n + 1), *y = t_jet(n + 1), *z = t_jet(n + 1);
    t_args(argv, argc, x, y, z, &b);
    long double *sx = t_jet(n), *sy = t_jet(n), *sz = t_jet(n), *cx = t_jet(n), *cy = t_jet(n), *cz = t_jet(n);

    // main loop
    t_xyz_output(x[0], y[0], z[0], 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            x[k + 1] = (*t_sin_cos(sy, cy, y, k, TRIG).a - b * x[k]) / (k + 1);
            y[k + 1] = (*t_sin_cos(sz, cz, z, k, TRIG).a - b * y[k]) / (k + 1);
            z[k + 1] = (*t_sin_cos(sx, cx, x, k, TRIG).a - b * z[k]) / (k + 1);
        }

        // sum the series using Horner's method and advance one step
        t_xyz_output(t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step);
    }
    return 0;
}
