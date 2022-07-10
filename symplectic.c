/*
 * Symplectic Integrators and support functions
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"

controls *get_c (char **argv) {
    controls *c = malloc(sizeof (controls));
    c->order = strtol(argv[2], NULL, 10); assert(c->order >= 2 && c->order <= 10);
    c->step_size = strtold(argv[3], NULL); assert(c->step_size > 0.0L);
    c->steps = strtol(argv[4], NULL, 10); assert(c->steps >= 0 && c->steps <= 1000000);
    return c;
}

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
} weight_r_4, weight_r_3, weight_r_2, weight_r_1;

static void stormer_verlet (void *p, real cd) {
    update_q(p, cd * 0.5L);
    update_p(p, cd);
    update_q(p, cd * 0.5L);
}

static void suzuki (void *p, integrator base, real cd, struct weight w) {
    base(p, cd * w.fwd);
    base(p, cd * w.fwd);
    base(p, cd * w.rev);
    base(p, cd * w.fwd);
    base(p, cd * w.fwd);
}

static void base4 (void *p, real cd) {
    suzuki(p, stormer_verlet, cd, weight_r_1);
}

static void fourth_order (void *p, real h) {
    base4(p, h);
}

static void base6 (void *p, real cd) {
    suzuki(p, base4, cd, weight_r_2);
}

static void sixth_order (void *p, real h) {
    base6(p, h);
}

static void base8 (void *p, real cd) {
    suzuki(p, base6, cd, weight_r_3);
}

static void eightth_order (void *p, real h) {
    base8(p, h);
}

static void base10 (void *p, real cd) {
    suzuki(p, base8, cd, weight_r_4);
}

static void tenth_order (void *p, real h) {
    base10(p, h);
}

static void set_weights (void) {
    weight_r_1.fwd = 1.0L / (4.0L - powl(4.0L, 1.0L / 3.0L));
    weight_r_2.fwd = 1.0L / (4.0L - powl(4.0L, 1.0L / 5.0L));
    weight_r_3.fwd = 1.0L / (4.0L - powl(4.0L, 1.0L / 7.0L));
    weight_r_4.fwd = 1.0L / (4.0L - powl(4.0L, 1.0L / 9.0L));
    weight_r_1.rev = 1.0L - 4.0L * weight_r_1.fwd;
    weight_r_2.rev = 1.0L - 4.0L * weight_r_2.fwd;
    weight_r_3.rev = 1.0L - 4.0L * weight_r_3.fwd;
    weight_r_4.rev = 1.0L - 4.0L * weight_r_4.fwd;
}

static integrator set_integrator (long order) {
    integrator composer = NULL;
    switch (order) {
        case 2: composer = stormer_verlet; break;
        case 4: composer = fourth_order; break;
        case 6: composer = sixth_order; break;
        case 8: composer = eightth_order; break;
        case 10: composer = tenth_order; break;
        default:
            printf("Order parameter is {%ld} but should be 2, 4, 6, 8, or 10 \n", order);
            exit(1);
    }
    return composer;
}

void solve (char **argv, void *p, plotter output) {
    int display_precision = (int)strtol(argv[1], NULL, 10); assert(display_precision >= 1 && display_precision <= 32);
    long order = strtol(argv[2], NULL, 10); assert(order >= 2 && order <= 10);
    real step_size = strtold(argv[3], NULL); assert(step_size > 0.0L);
    long steps = strtol(argv[4], NULL, 10); assert(steps >= 0 && steps <= 1000000);
    set_weights();
    integrator composer = set_integrator(order);
    for (long step = 0; step < steps; step++) {
        output(display_precision, p, step * step_size);
        composer(p, step_size);
    }
    output(display_precision, p, steps * step_size);
}

void *generate (controls *cont, void *params) {
    static controls *c;
    static void *p;
    static integrator composer = NULL;
    static long step, resume = 0;
    if (resume) goto resume; else resume = 1;
    c = cont;
    p = params;
    set_weights();
    composer = set_integrator(c->order);
    for (step = 1; step <= c->steps; step++) {
        composer(p, c->step_size);
        return p;
        resume: ;
    }
    resume = 0;
    return NULL;
}
