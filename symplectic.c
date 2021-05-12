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

void t_stepper (char **argv, long *dp, long *method, real *h, long *nsteps) {
    *dp = strtol(argv[1], NULL, 10);
    *method = strtol(argv[2], NULL, 10);
    *h = strtold(argv[3], NULL);
    assert(*h > 0.0L && *h <= 10.0L);
    *nsteps = strtol(argv[4], NULL, 10);
    assert(*nsteps >= 0 && *nsteps <= 1000000);
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

static struct {
    real w0, x0, y0, z0, w1, x1, y1, z1;
} weight;

static void stormer_verlet (void *p, updater uq, updater up, real cd) {
    uq(p, cd * 0.5L);
    up(p, cd);
    uq(p, cd * 0.5L);
}

static void suzuki (void *p, integrator base, updater uq, updater up, real cd, real forward, real back) {
    base(p, uq, up, cd * forward);
    base(p, uq, up, cd * forward);
    base(p, uq, up, cd * back);
    base(p, uq, up, cd * forward);
    base(p, uq, up, cd * forward);
}

static void base4_suzuki (void *p, updater uq, updater up, real cd) {
    suzuki(p, stormer_verlet, uq, up, cd, weight.z1, weight.z0);
}

static void fourth_order_suzuki (void *p, updater uq, updater up, real h) {
    base4_suzuki(p, uq, up, h);
}

static void base6_suzuki (void *p, updater uq, updater up, real cd) {
    suzuki(p, base4_suzuki, uq, up, cd, weight.y1, weight.y0);
}

static void sixth_order_suzuki (void *p, updater uq, updater up, real h) {
    base6_suzuki(p, uq, up, h);
}

static void base8_suzuki (void *p, updater uq, updater up, real cd) {
    suzuki(p, base6_suzuki, uq, up, cd, weight.x1, weight.x0);
}

static void eightth_order_suzuki (void *p, updater uq, updater up, real h) {
    base8_suzuki(p, uq, up, h);
}

static void base10_suzuki (void *p, updater uq, updater up, real cd) {
    suzuki(p, base8_suzuki, uq, up, cd, weight.w1, weight.w0);
}

static void tenth_order_suzuki (void *p, updater uq, updater up, real h) {
    base10_suzuki(p, uq, up, h);
}

void solve (char **argv, void *p, updater uq, updater up, plotter output) {
    long method, steps, dp;
    real h;
    t_stepper(argv, &dp, &method, &h, &steps);
    weight.z1 = 1.0L / (4.0L - powl(4.0L, (1.0L / 3.0L)));
    weight.y1 = 1.0L / (4.0L - powl(4.0L, (1.0L / 5.0L)));
    weight.x1 = 1.0L / (4.0L - powl(4.0L, (1.0L / 7.0L)));
    weight.w1 = 1.0L / (4.0L - powl(4.0L, (1.0L / 9.0L)));
    weight.z0 = 1.0L - 4.0L * weight.z1;
    weight.y0 = 1.0L - 4.0L * weight.y1;
    weight.x0 = 1.0L - 4.0L * weight.x1;
    weight.w0 = 1.0L - 4.0L * weight.w1;
    integrator composer = NULL;
    switch (method) {
        case 2: composer = stormer_verlet; break;
        case 4: composer = fourth_order_suzuki; break;
        case 6: composer = sixth_order_suzuki; break;
        case 8: composer = eightth_order_suzuki; break;
        case 10: composer = tenth_order_suzuki; break;
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
