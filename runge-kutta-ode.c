/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "runge-kutta-ode.h"

static char *tag (real new, real old, char *max) {
    return new * old < 0.0L ? max : "_";
}

void rk4 (int dp, int n, real h, int steps, real x0, real y0, real z0, void *p, clock_t t0) {
    real x = x0, x_1 = 0.0L;
    real y = y0, y_1 = 0.0L;
    real z = z0, z_1 = 0.0L;
    components s = (components) {0.0L, 0.0L, 0.0L};
    for (int step = 0; step < steps; step++) {
        if (step % n == 0) {
            t_out(dp, x, y, z, h * step, tag(x_1, s.x, "X"), tag(y_1, s.y, "Y"), tag(z_1, s.z, "Z"), t0);
            s = (components) {x_1, y_1, z_1};
        }
        components k1 = ode(x, y, z, p);
        components k2 = ode(x + 0.5L * k1.x * h, y + 0.5L * k1.y * h, z + 0.5L * k1.z * h, p);
        components k3 = ode(x + 0.5L * k2.x * h, y + 0.5L * k2.y * h, z + 0.5L * k2.z * h, p);
        components k4 = ode(x + k3.x * h, y + k3.y * h, z + k3.z * h, p);
        x += h * (x_1 = (k1.x + 2.0L * (k2.x + k3.x) + k4.x) / 6.0L);
        y += h * (y_1 = (k1.y + 2.0L * (k2.y + k3.y) + k4.y) / 6.0L);
        z += h * (z_1 = (k1.z + 2.0L * (k2.z + k3.z) + k4.z) / 6.0L);
        if (! isfinite(x) || ! isfinite(y) || ! isfinite(z)) {
            fprintf(stderr, "Value error!\n");
            exit(2);
        }
    }
    t_out(dp, x, y, z, h * steps, "_", "_", "_", t0);
}
