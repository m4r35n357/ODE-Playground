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

void rk4 (int dp, controls *c, components *coordinates, void *p, clock_t t0) {
    real x = coordinates->x, x_1 = 0.0L;
    real y = coordinates->y, y_1 = 0.0L;
    real z = coordinates->z, z_1 = 0.0L;
    components s = (components) {0.0L, 0.0L, 0.0L};
    for (int step = 0; step < c->steps; step++) {
        if (step % c->order == 0) {
            t_out(dp, x, y, z, c->step_size * step, tag(x_1, s.x, "X"), tag(y_1, s.y, "Y"), tag(z_1, s.z, "Z"), t0);
            s = (components) {x_1, y_1, z_1};
        }
        components k1 = ode(x, y, z, p);
        components k2 = ode(x + 0.5L * k1.x * c->step_size, y + 0.5L * k1.y * c->step_size, z + 0.5L * k1.z * c->step_size, p);
        components k3 = ode(x + 0.5L * k2.x * c->step_size, y + 0.5L * k2.y * c->step_size, z + 0.5L * k2.z * c->step_size, p);
        components k4 = ode(x + k3.x * c->step_size, y + k3.y * c->step_size, z + k3.z * c->step_size, p);
        x += c->step_size * (x_1 = (k1.x + 2.0L * (k2.x + k3.x) + k4.x) / 6.0L);
        y += c->step_size * (y_1 = (k1.y + 2.0L * (k2.y + k3.y) + k4.y) / 6.0L);
        z += c->step_size * (z_1 = (k1.z + 2.0L * (k2.z + k3.z) + k4.z) / 6.0L);
        if (! isfinite(x) || ! isfinite(y) || ! isfinite(z)) {
            fprintf(stderr, "Value error!\n");
            exit(2);
        }
    }
    t_out(dp, x, y, z, c->step_size * c->steps, "_", "_", "_", t0);
}
