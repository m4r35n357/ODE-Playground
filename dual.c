
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "dual.h"

dual d_dual (real a) {
    return (dual) { .val = a, .dot = 0.0L };
}

dual d_var (real a) {
    return (dual) { .val = a, .dot = 1.0L };
}

dual d_abs (dual a) {
    return (dual) { .val = a.val < 0.0L ? - a.val : a.val, .dot = a.val < 0.0L ? - a.dot : a.dot };
}

dual d_neg (dual b) {
    return (dual) { .val = - b.val, .dot = - b.dot };
}

dual d_inv (dual b) {
    assert(b.val != 0.0L);
    return (dual) { .val = 1.0L / b.val, .dot = - b.dot / (b.val * b.val) };
}

dual d_sqr (dual a) {
    return (dual) { .val = a.val * a.val, .dot = 2.0L * a.val * a.dot };
}

dual d_shift (dual a, real b) {
    return (dual) { .val = a.val + b, .dot = a.dot };
}

dual d_scale (dual a, real b) {
    return (dual) { .val = a.val * b, .dot = a.dot * b };
}

dual d_add (dual a, dual b) {
    return (dual) { .val = a.val + b.val, .dot = a.dot + b.dot };
}

dual d_sub (dual a, dual b) {
    return (dual) { .val = a.val - b.val, .dot = a.dot - b.dot };
}

dual d_mul (dual a, dual b) {
    return (dual) { .val = a.val * b.val, .dot = a.val * b.dot + a.dot * b.val };
}

dual d_div (dual a, dual b) {
    assert(b.val != 0.0L);
    return (dual) { .val = a.val / b.val, .dot = (a.dot * b.val - a.val * b.dot) / (b.val * b.val) };
}

dual d_exp (dual a) {
    real exp_val = expl(a.val);
    return (dual) { .val = exp_val, .dot = a.dot * exp_val };
}

dual d_log (dual a) {
    assert(a.val > 0.0L);
    return (dual) { .val = logl(a.val), .dot = a.dot / a.val };
}

dual d_pow (dual a, real b) {
    assert(a.val > 0.0L);
    return (dual) { .val = powl(a.val, b), .dot = a.dot * b * powl(a.val, (b - 1)) };
}

dual d_sin (dual a) {
    return (dual) { .val = sinl(a.val), .dot = a.dot * cosl(a.val) };
}

dual d_cos (dual a) {
    return (dual) { .val = cosl(a.val), .dot = - a.dot * sinl(a.val) };
}

dual d_tan (dual a) {
    real tan_val = tanl(a.val);
    return (dual) { .val = tan_val, .dot = a.dot * (1.0L + tan_val * tan_val) };
}

dual d_sinh (dual a) {
    return (dual) { .val = sinhl(a.val), .dot = a.dot * coshl(a.val) };
}

dual d_cosh (dual a) {
    return (dual) { .val = coshl(a.val), .dot = a.dot * sinhl(a.val) };
}

dual d_tanh (dual a) {
    real tanh_val = tanhl(a.val);
    return (dual) { .val = tanh_val, .dot = a.dot * (1.0L - tanh_val * tanh_val) };
}

void t_stepper (char **argv, long *dp, long *method, real *h, long *nsteps) {
    *dp = strtol(argv[1], NULL, 10);
    *method = strtol(argv[2], NULL, 10);
    *h = strtold(argv[3], NULL);
    assert(*h > 0.0L && *h <= 10.0L);
    *nsteps = strtol(argv[4], NULL, 10);
    assert(*nsteps >= 1 && *nsteps <= 100000);
}

void t_variables (char **argv, int argc, ...) {
    va_list vars;
    va_start(vars, argc);
    for (int i = 5; i < argc; i++) {
        *va_arg(vars, real *) = strtold(argv[i], NULL);
    }
    va_end(vars);
}

static void cromer (void *p, long size, real cd[], updater uq, updater up) {
    (void)size;
    uq(p, cd[0]);
    up(p, cd[1]);
}

static void symplectic (void *p, long size, real cd[], updater uq, updater up) {
    for (long i = 0; i < size; i++) {
        i % 2 == 0 ? uq(p, cd[i]) : up(p, cd[i]);  // 0 -> size -1
    }
    for (long i = size - 2; i > -1; i--) {
        i % 2 == 0 ? uq(p, cd[i]) : up(p, cd[i]);  // size -2 -> 0
    }
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
    integrator i;
    switch (method) {
        case 1:  // Cromer's method (symplectic Euler)
            i = cromer;
            size = 2;
            cd = (real[]){ h, h };
            break;
        case 2:  // Stormer-Verlet
            i = symplectic;
            size = 2;
            cd = (real[]){ 0.5L * h, h };
            break;
        case -4:  // Forest-Ruth (Yoshida-style composition)
            i = symplectic;
            size = 4;
            z1 = 1.0L / (2.0L - powl(2.0L, (1.0L / 3.0L)));
            z0 = 1.0L - 2.0L * z1;
            cd = (real[]){0.5L * h * z1, h * z1, 0.5L * h * (z1 + z0), h * z0};
            break;
        case 4:  // Like Forest-Ruth but with Suzuki-style composition
            i = symplectic;
            size = 6;
            cd = (real[]){0.5L * h * z1, h * z1, h * z1, h * z1, 0.5L * h * (z1 + z0), h * z0};
            break;
        case 6:  // Higher order Suzuki
            i = symplectic;
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
                0.5L * h * z1 * (y1 + y0),
                h * z1 * y0, h * z1 * y0, h * z1 * y0,
                0.5L * h * (z1 + z0) * y0, h * z0 * y0
            };
            break;
        case 8:  // Higher order Suzuki
            i = symplectic;
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
            printf("Method parameter is {%ld} but should be 1, 2, 4, 6, 8 or -4\n", method);
            exit(1);
    }
    for (long step = 1; step <= steps; step++) {
        i(p, size, cd, uq, up);
        output(dp, p, step * h);
    }
}
