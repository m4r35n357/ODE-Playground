/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

void tsm (int argc, char **argv, long dp, long n, real h, long steps, real x0, real y0, real z0) {
    series x = t_jet(n + 1), y = t_jet(n + 1), z = t_jet(n + 1);
    x[0] = x0; y[0] = y0; z[0] = z0;
    void *p = get_p(argc, argv, n);
    components cdot = ode(x, y, z, p, 0);
    t_output(dp, x[0], y[0], z[0], 0.0L, "_", "_", "_");
    for (long step = 1; step <= steps; step++) {
        for (int k = 0; k < n; k++) {
            components c = ode(x, y, z, p, k);
            x[k + 1] = c.x / (k + 1);
            y[k + 1] = c.y / (k + 1);
            z[k + 1] = c.z / (k + 1);
        }
        t_output(dp, t_horner(x, n, h), t_horner(y, n, h), t_horner(z, n, h), h * step,
                 x[1] * cdot.x < 0.0L ? (x[2] > 0.0L ? "x" : "X") : "_",
                 y[1] * cdot.y < 0.0L ? (y[2] > 0.0L ? "y" : "Y") : "_",
                 z[1] * cdot.z < 0.0L ? (z[2] > 0.0L ? "z" : "Z") : "_");
        cdot.x = x[1]; cdot.y = y[1]; cdot.z = z[1];
    }
}

int main (int argc, char **argv) {
    long dp = strtol(argv[1], NULL, 10); assert(dp >= 1 && dp <= 99);
    long n = strtol(argv[2], NULL, 10); assert(n >= 2 && n <= 64);
    real h = strtold(argv[3], NULL); assert(h > 0.0L);
    long steps = strtol(argv[4], NULL, 10); assert(steps >= 1 && steps <= 1000000);

    tsm(argc, argv, dp, n, h, steps, strtold(argv[5], NULL), strtold(argv[6], NULL), strtold(argv[7], NULL));

    return 0;
}
