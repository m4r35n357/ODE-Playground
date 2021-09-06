/*
 * Symplectic Integrators and support functions
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"

void t_variables (char **argv, int begin, int argc, ...) {
    va_list model_params;
    va_start(model_params, argc);
    for (int i = begin; i < argc; i++) {
        *va_arg(model_params, real *) = strtold(argv[i], NULL);
    }
    va_end(model_params);
}

real error (real e) {
    return 10.0L * log10l(fabsl(e) >= 1e-36L ? fabsl(e) : 1e-36L);
}

static struct weight {
    real fwd, rev;
} w_w, w_x, w_y, w_z;

static void stormer_verlet (void *p, updater uq, updater up, real cd) {
    uq(p, cd * 0.5L);
    up(p, cd);
    uq(p, cd * 0.5L);
}

static void suzuki (void *p, integrator base, updater uq, updater up, real cd, struct weight w) {
    base(p, uq, up, cd * w.fwd);
    base(p, uq, up, cd * w.fwd);
    base(p, uq, up, cd * w.rev);
    base(p, uq, up, cd * w.fwd);
    base(p, uq, up, cd * w.fwd);
}

static void base4 (void *p, updater uq, updater up, real cd) {
    suzuki(p, stormer_verlet, uq, up, cd, w_z);
}

static void fourth_order (void *p, updater uq, updater up, real h) {
    base4(p, uq, up, h);
}

static void base6 (void *p, updater uq, updater up, real cd) {
    suzuki(p, base4, uq, up, cd, w_y);
}

static void sixth_order (void *p, updater uq, updater up, real h) {
    base6(p, uq, up, h);
}

static void base8 (void *p, updater uq, updater up, real cd) {
    suzuki(p, base6, uq, up, cd, w_x);
}

static void eightth_order (void *p, updater uq, updater up, real h) {
    base8(p, uq, up, h);
}

static void base10 (void *p, updater uq, updater up, real cd) {
    suzuki(p, base8, uq, up, cd, w_w);
}

static void tenth_order (void *p, updater uq, updater up, real h) {
    base10(p, uq, up, h);
}

void solve (char **argv, void *p, updater uq, updater up, plotter output) {
    long dp = strtol(argv[1], NULL, 10), method = strtol(argv[2], NULL, 10), steps = strtol(argv[4], NULL, 10);
    real h = strtold(argv[3], NULL);
    assert(h > 0.0L && h <= 10.0L);
    assert(steps >= 0 && steps <= 1000000);
    w_z.fwd = 1.0L / (4.0L - powl(4.0L, 1.0L / 3.0L));
    w_y.fwd = 1.0L / (4.0L - powl(4.0L, 1.0L / 5.0L));
    w_x.fwd = 1.0L / (4.0L - powl(4.0L, 1.0L / 7.0L));
    w_w.fwd = 1.0L / (4.0L - powl(4.0L, 1.0L / 9.0L));
    w_z.rev = 1.0L - 4.0L * w_z.fwd;
    w_y.rev = 1.0L - 4.0L * w_y.fwd;
    w_x.rev = 1.0L - 4.0L * w_x.fwd;
    w_w.rev = 1.0L - 4.0L * w_w.fwd;
    integrator composer = NULL;
    switch (method) {
        case 2: composer = stormer_verlet; break;
        case 4: composer = fourth_order; break;
        case 6: composer = sixth_order; break;
        case 8: composer = eightth_order; break;
        case 10: composer = tenth_order; break;
        default:
            printf("Method parameter is {%ld} but should be 2, 4, 6, 8, or 10 \n", method);
            exit(1);
    }
    output(dp, p, 0.0L);
    for (long step = 1; step <= steps; step++) {
        composer(p, uq, up, h);
        output(dp, p, step * h);
    }
}
