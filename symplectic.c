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

struct stage {
    real w0, x0, y0, z0, w1, x1, y1, z1;
} s;

static void euler_cromer (void *p, updater uq, updater up, real cd) {
    uq(p, cd);
    up(p, cd);
}

static void base2 (void *p, updater uq, updater up, real cd) {
    uq(p, cd * 0.5L);
    up(p, cd);
    uq(p, cd * 0.5L);
}

static void second_order (void *p, updater uq, updater up, real h) {
    base2(p, uq, up, h);
}

static void yoshida (void *p, base b, updater uq, updater up, real cd, real forward, real back) {
    b(p, uq, up, cd * forward);
    b(p, uq, up, cd * back);
    b(p, uq, up, cd * forward);
}

static void base4_yoshida (void *p, updater uq, updater up, real cd) {
    yoshida(p, base2, uq, up, cd, s.z1, s.z0);
}

static void fourth_order_yoshida (void *p, updater uq, updater up, real h) {
    base4_yoshida(p, uq, up, h);
}

static void base6_yoshida (void *p, updater uq, updater up, real cd) {
    yoshida(p, base4_yoshida, uq, up, cd, s.y1, s.y0);
}

static void sixth_order_yoshida (void *p, updater uq, updater up, real h) {
    base6_yoshida(p, uq, up, h);
}

static void base8_yoshida (void *p, updater uq, updater up, real cd) {
    yoshida(p, base6_yoshida, uq, up, cd, s.x1, s.x0);
}

static void eightth_order_yoshida (void *p, updater uq, updater up, real h) {
    base8_yoshida(p, uq, up, h);
}

static void base10_yoshida (void *p, updater uq, updater up, real cd) {
    yoshida(p, base8_yoshida, uq, up, cd, s.w1, s.w0);
}

static void tenth_order_yoshida (void *p, updater uq, updater up, real h) {
    base10_yoshida(p, uq, up, h);
}

static void suzuki (void *p, base b, updater uq, updater up, real cd, real forward, real back) {
    b(p, uq, up, cd * forward);
    b(p, uq, up, cd * forward);
    b(p, uq, up, cd * back);
    b(p, uq, up, cd * forward);
    b(p, uq, up, cd * forward);
}

static void base4_suzuki (void *p, updater uq, updater up, real cd) {
    suzuki(p, base2, uq, up, cd, s.z1, s.z0);
}

static void fourth_order_suzuki (void *p, updater uq, updater up, real h) {
    base4_suzuki(p, uq, up, h);
}

static void base6_suzuki (void *p, updater uq, updater up, real cd) {
    suzuki(p, base4_suzuki, uq, up, cd, s.y1, s.y0);
}

static void sixth_order_suzuki (void *p, updater uq, updater up, real h) {
    base6_suzuki(p, uq, up, h);
}

static void base8_suzuki (void *p, updater uq, updater up, real cd) {
    suzuki(p, base6_suzuki, uq, up, cd, s.x1, s.x0);
}

static void eightth_order_suzuki (void *p, updater uq, updater up, real h) {
    base8_suzuki(p, uq, up, h);
}

static void base10_suzuki (void *p, updater uq, updater up, real cd) {
    suzuki(p, base8_suzuki, uq, up, cd, s.w1, s.w0);
}

static void tenth_order_suzuki (void *p, updater uq, updater up, real h) {
    base10_suzuki(p, uq, up, h);
}

void solve (char **argv, void *p, updater uq, updater up, plotter output) {
    long method, steps, dp;
    real h;
    t_stepper(argv, &dp, &method, &h, &steps);
    integrator composer = NULL;
    if (method < 0L) {
        s.z1 = 1.0L / (2.0L - powl(2.0L, (1.0L / 3.0L)));
        s.y1 = 1.0L / (2.0L - powl(2.0L, (1.0L / 5.0L)));
        s.x1 = 1.0L / (2.0L - powl(2.0L, (1.0L / 7.0L)));
        s.w1 = 1.0L / (2.0L - powl(2.0L, (1.0L / 9.0L)));
        s.z0 = 1.0L - 2.0L * s.z1;
        s.y0 = 1.0L - 2.0L * s.y1;
        s.x0 = 1.0L - 2.0L * s.x1;
        s.w0 = 1.0L - 2.0L * s.w1;
        switch (method) {
            case -4:  // Composed Forest-Ruth
                composer = fourth_order_yoshida;
                break;
            case -6:  // Higher order Yoshida
                composer = sixth_order_yoshida;
                break;
            case -8:  // Higher order Yoshida
                composer = eightth_order_yoshida;
                break;
            case -10:  // Higher order Yoshida
                composer = tenth_order_yoshida;
                break;
            default:
                printf("Method parameter is {%ld} but should be -4, -6, -8, or -10\n", method);
                exit(1);
        }
    } else {
        s.z1 = 1.0L / (4.0L - powl(4.0L, (1.0L / 3.0L)));
        s.y1 = 1.0L / (4.0L - powl(4.0L, (1.0L / 5.0L)));
        s.x1 = 1.0L / (4.0L - powl(4.0L, (1.0L / 7.0L)));
        s.w1 = 1.0L / (4.0L - powl(4.0L, (1.0L / 9.0L)));
        s.z0 = 1.0L - 4.0L * s.z1;
        s.y0 = 1.0L - 4.0L * s.y1;
        s.x0 = 1.0L - 4.0L * s.x1;
        s.w0 = 1.0L - 4.0L * s.w1;
        switch (method) {
            case 1:  // Euler-Cromer
                composer = euler_cromer;
                break;
            case 2:  // Stormer-Verlet
                composer = second_order;
                break;
            case 4:  // Like composed Forest-Ruth but with Suzuki-style composition
                composer = fourth_order_suzuki;
                break;
            case 6:  // Higher order Suzuki
                composer = sixth_order_suzuki;
                break;
            case 8:  // Higher order Suzuki
                composer = eightth_order_suzuki;
                break;
            case 10:  // Higher order Suzuki
                composer = tenth_order_suzuki;
                break;
            default:
                printf("Method parameter is {%ld} but should be 1, 2, 4, 6, 8, or 10 \n", method);
                exit(1);
        }
    }
    output(dp, p, 0.0L);
    for (long step = 1; step <= steps; step++) {
        composer(p, uq, up, h);
        output(dp, p, step * h);
    }
}
