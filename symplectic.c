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

static real r (int stage) {
    return 1.0L / (4.0L - powl(4.0L, 1.0L / (2.0L * (real)stage + 1.0L)));
}

controls *get_c (char **argv) {
    controls *c = malloc(sizeof (controls));
    c->order = strtol(argv[2], NULL, 10); assert(c->order >= 2 && c->order <= 10);
    c->step_size = strtold(argv[3], NULL); assert(c->step_size > 0.0L);
    c->steps = strtol(argv[4], NULL, 10); assert(c->steps >= 0 && c->steps <= 1000000);
    c->r1 = (weights) { .fwd = r(1), .rev = 1.0L - 4.0L * r(1) };
    c->r2 = (weights) { .fwd = r(2), .rev = 1.0L - 4.0L * r(2) };
    c->r3 = (weights) { .fwd = r(3), .rev = 1.0L - 4.0L * r(3) };
    c->r4 = (weights) { .fwd = r(4), .rev = 1.0L - 4.0L * r(4) };
    return c;
}

real error (real e) {
    return 10.0L * log10l(fabsl(e) >= 1e-36L ? fabsl(e) : 1e-36L);
}

static void stormer_verlet (controls *cont, void *p, real cd) { (void)cont;
    update_q(p, cd * 0.5L);
    update_p(p, cd);
    update_q(p, cd * 0.5L);
}

static void suzuki (controls *c, void *p, integrator base, real cd, weights w) {
    base(c, p, cd * w.fwd);
    base(c, p, cd * w.fwd);
    base(c, p, cd * w.rev);
    base(c, p, cd * w.fwd);
    base(c, p, cd * w.fwd);
}

static void base4 (controls *c, void *p, real cd) {
    suzuki(c, p, stormer_verlet, cd, c->r1);
}

static void fourth_order (controls *c, void *p, real h) {
    base4(c, p, h);
}

static void base6 (controls *c, void *p, real cd) {
    suzuki(c, p, base4, cd, c->r2);
}

static void sixth_order (controls *c, void *p, real h) {
    base6(c, p, h);
}

static void base8 (controls *c, void *p, real cd) {
    suzuki(c, p, base6, cd, c->r3);
}

static void eightth_order (controls *c, void *p, real h) {
    base8(c, p, h);
}

static void base10 (controls *c, void *p, real cd) {
    suzuki(c, p, base8, cd, c->r4);
}

static void tenth_order (controls *c, void *p, real h) {
    base10(c, p, h);
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
    controls *c = get_c(argv);
    integrator composer = set_integrator(c->order);
    for (long step = 0; step < c->steps; step++) {
        output(display_precision, p, step * c->step_size);
        composer(c, p, c->step_size);
    }
    output(display_precision, p, c->steps * c->step_size);
}

void *generate (controls *cont, void *params) {
    static controls *c;
    static void *p;
    static integrator composer = NULL;
    static long step, resume = 0;
    if (resume) goto resume; else resume = 1;
    c = cont;
    p = params;
    composer = set_integrator(c->order);
    for (step = 1; step <= c->steps; step++) {
        c->step = step;
        composer(c, p, c->step_size);
        return p;
        resume: ;
    }
    resume = 0;
    return NULL;
}
