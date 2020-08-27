/*
 * Sprott-Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-sprott-thomas-dbg NA NA 10 0.01 30000 1 0 0 4.75 1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

int main (int argc, char **argv) {
    long n, nsteps;
    long double a, b, h;

    // initialize from command arguments
    assert(argc == 11);
    t_stepper(argv, &n, &h, &nsteps);
    long double *x = t_jet(n + 1), *y = t_jet(n + 1), *z = t_jet(n + 1);
    t_args(argv, argc, x, y, z, &a, &b);
    long double *ax = t_jet(n), *ay = t_jet(n), *az = t_jet(n);
    long double *sax = t_jet(n), *say = t_jet(n), *saz = t_jet(n), *cax = t_jet(n), *cay = t_jet(n), *caz = t_jet(n);
    long double *tx = t_jet(n), *ty = t_jet(n), *tz = t_jet(n), *s2x = t_jet(n), *s2y = t_jet(n), *s2z = t_jet(n);

    // main loop
    t_xyz_output(x[0], y[0], z[0], 0.0);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            ax[k] = a * x[k];
            ay[k] = a * y[k];
            az[k] = a * z[k];
            x[k + 1] = (*t_sin_cos(say, cay, ay, k, TRIG).a - b * *t_tan_sec2(tx, s2x, x, k, TRIG).a) / (k + 1);
            y[k + 1] = (*t_sin_cos(saz, caz, az, k, TRIG).a - b * *t_tan_sec2(ty, s2y, y, k, TRIG).a) / (k + 1);
            z[k + 1] = (*t_sin_cos(sax, cax, ax, k, TRIG).a - b * *t_tan_sec2(tz, s2z, z, k, TRIG).a) / (k + 1);
        }

        // sum the series using Horner's method and advance one step
        t_xyz_output(t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step);
    }
    return 0;
}
