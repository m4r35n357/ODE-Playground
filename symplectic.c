
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"

void t_stepper (char **argv, long *dp, long *method, real *h, long *nsteps) {
    *dp = strtol(argv[1], NULL, 10);
    *method = strtol(argv[2], NULL, 10);
    *h = strtold(argv[3], NULL);
    assert(*h > 0.0L && *h <= 10.0L);
    *nsteps = strtol(argv[4], NULL, 10);
    assert(*nsteps >= 1 && *nsteps <= 1000000);
}

void t_variables (char **argv, int begin, int argc, ...) {
    va_list vars;
    va_start(vars, argc);
    for (int i = begin; i < argc; i++) {
        *va_arg(vars, real *) = strtold(argv[i], NULL);
    }
    va_end(vars);
}

real error (real e) {
    return 10.0L * log10l(fabsl(e) >= 1e-36L ? fabsl(e) : 1e-36L);
}

void solve (char **argv, void *p, updater uq, updater up, plotter output) {
    long method, steps, dp, size;
    real h, x0, x1, y0, y1, z0, z1, *cd = (real[]){ 0.0L };
    t_stepper(argv, &dp, &method, &h, &steps);
    z1 = 1.0L / (4.0L - powl(4.0L, (1.0L / 3.0L)));
    y1 = 1.0L / (4.0L - powl(4.0L, (1.0L / 5.0L)));
    x1 = 1.0L / (4.0L - powl(4.0L, (1.0L / 7.0L)));
    z0 = 1.0L - 4.0L * z1;
    y0 = 1.0L - 4.0L * y1;
    x0 = 1.0L - 4.0L * x1;
    switch (method) {
        case 2:  // Stormer-Verlet
            size = 2;
            cd = (real[]){ 0.5L * h, h };
            break;
        case 4:  // Like Forest-Ruth but with Suzuki-style composition
            size = 6;
            cd = (real[]){0.5L * h * z1, h * z1, h * z1, h * z1, 0.5L * h * (z1 + z0), h * z0};
            break;
        case 6:  // Higher order Suzuki
            size = 26;
            cd = (real[]){
                0.5L * h * z1 * y1,
                h * z1 * y1, h * z1 * y1, h * z1 * y1,
                0.5L * h * (z1 + z0) * y1, h * z0 * y1, 0.5L * h * (z0 + z1) * y1,
                h * z1 * y1, h * z1 * y1, h * z1 * y1,
                h * z1 * y1,
                h * z1 * y1, h * z1 * y1, h * z1 * y1,
                0.5L * h * (z1 + z0) * y1, h * z0 * y1, 0.5L * h * (z0 + z1) * y1,
                h * z1 * y1, h * z1 * y1, h * z1 * y1,
                0.5L * h * z1 * (y1 + y0), h * z1 * y0, h * z1 * y0, h * z1 * y0, 0.5L * h * (z1 + z0) * y0, h * z0 * y0
            };
            break;
        case 8:  // Higher order Suzuki
            size = 126;
            cd = (real[]){
                0.5L * h * z1 * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                0.5L * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5L * h * (z0 + z1) * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                h * z1 * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                0.5L * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5L * h * (z0 + z1) * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                0.5L * h * z1 * (y1 + y0) * x1,
                h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                0.5L * h * (z1 + z0) * y0 * x1, h * z0 * y0 * x1, 0.5L * h * (z0 + z1) * y0 * x1,
                h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                0.5L * h * z1 * (y1 + y0) * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                0.5L * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5L * h * (z0 + z1) * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                h * z1 * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                0.5L * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5L * h * (z0 + z1) * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                h * z1 * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                0.5L * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5L * h * (z0 + z1) * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                h * z1 * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                0.5L * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5L * h * (z0 + z1) * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                0.5L * h * z1 * (y1 + y0) * x1,
                h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                0.5L * h * (z1 + z0) * y0 * x1, h * z0 * y0 * x1, 0.5L * h * (z0 + z1) * y0 * x1,
                h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                0.5L * h * z1 * (y1 + y0) * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                0.5L * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5L * h * (z0 + z1) * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                h * z1 * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                0.5L * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5L * h * (z0 + z1) * y1 * x1,
                h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                0.5L * h * z1 * y1 * (x1 + x0),
                h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                0.5L * h * (z1 + z0) * y1 * x0, h * z0 * y1 * x0, 0.5L * h * (z0 + z1) * y1 * x0,
                h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                h * z1 * y1 * x0,
                h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                0.5L * h * (z1 + z0) * y1 * x0, h * z0 * y1 * x0, 0.5L * h * (z0 + z1) * y1 * x0,
                h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                0.5L * h * z1 * (y1 + y0) * x0,
                h * z1 * y0 * x0, h * z1 * y0 * x0, h * z1 * y0 * x0,
                0.5L * h * (z1 + z0) * y0 * x0, h * z0 * y0 * x0
            };
            break;
        default:
            printf("Method parameter is {%ld} but should be 2, 4, 6, or 8\n", method);
            exit(1);
    }
    output(dp, p, 0.0L);
    for (long step = 1; step <= steps; step++) {
        for (long i = 0; i < size; i++) {
            i % 2 == 0 ? uq(p, cd[i]) : up(p, cd[i]);  // 0 -> size -1
        }
        for (long i = size - 2; i > -1; i--) {
            i % 2 == 0 ? uq(p, cd[i]) : up(p, cd[i]);  // size -2 -> 0
        }
        output(dp, p, step * h);
    }
}
