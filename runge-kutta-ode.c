/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "runge-kutta-ode.h"

static char *tag (real new, real old, char *max) {
    return new * old < 0.0L ? max : "_";
}

void rk4 (int dp, int n, real h, int steps, real x0, real y0, real z0, void *p) {
    real x = x0, xdot = 0.0L, kx = 0.0L;
    real y = y0, ydot = 0.0L, ky = 0.0L;
    real z = z0, zdot = 0.0L, kz = 0.0L;
    clock_t start = clock();
    for (long step = 1; step < steps + 1; step++) {
        components k1 = ode(x, y, z, p);
        components k2 = ode(x + 0.5L * k1.x * h, y + 0.5L * k1.y * h, z + 0.5L * k1.z * h, p);
        components k3 = ode(x + 0.5L * k2.x * h, y + 0.5L * k2.y * h, z + 0.5L * k2.z * h, p);
        components k4 = ode(x + k3.x * h, y + k3.y * h, z + k3.z * h, p);
        x += h * (kx = (k1.x + 2.0L * (k2.x + k3.x) + k4.x) / 6.0L);
        y += h * (ky = (k1.y + 2.0L * (k2.y + k3.y) + k4.y) / 6.0L);
        z += h * (kz = (k1.z + 2.0L * (k2.z + k3.z) + k4.z) / 6.0L);
        if (! isfinite(x) || ! isfinite(y) || ! isfinite(z)) {
            fprintf(stderr, "Value error!\n");
            exit(2);
        }
        if (step % n == 0) {
            t_output(dp, x, y, z, h * step, tag(kx, xdot, "X"), tag(ky, ydot, "Y"), tag(kz, zdot, "Z"), cpu(start));
            xdot = kx, ydot = ky, zdot = kz;
        }
    }
    t_output(dp, x, y, z, h * steps, "_", "_", "_", cpu(start));
}
