/*
 * Symplectic Integrators and support functions
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symplectic.h"

const int BASE = 10;

controls *get_c_symp (int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    controls *c = malloc(sizeof (controls)); CHECK(c);
    c->order = (int)strtol(argv[2], NULL, BASE); CHECK(c->order >= 2 && c->order <= 10 && c->order % 2 == 0);
    c->step_size = strtold(argv[3], NULL);       CHECK(c->step_size > 0.0L);
    c->steps = (int)strtol(argv[4], NULL, BASE); CHECK(c->steps >= 0 && c->steps <= 1000000);
    return c;
}

real error (real e) {
    return - log10l(fabsl(e) >= 1e-36L ? fabsl(e) : 1e-36L);
}

static void second_order (controls *cont, void *p, real cd) { (void)cont;
    update_q(p, cd * 0.5L);
    update_p(p, cd);
    update_q(p, cd * 0.5L);
}

static void suzuki (controls *c, void *p, integrator base, real cd, weights w) {
    for (int step = 0; step < 5; step++) base(c, p, cd * (step == 2 ? w.rev : w.fwd));
}

static void base4 (controls *c, void *p, real cd) {
    suzuki(c, p, second_order, cd, c->r2);
}

static void fourth_order (controls *c, void *p, real h) {
    base4(c, p, h);
}

static void base6 (controls *c, void *p, real cd) {
    suzuki(c, p, base4, cd, c->r4);
}

static void sixth_order (controls *c, void *p, real h) {
    base6(c, p, h);
}

static void base8 (controls *c, void *p, real cd) {
    suzuki(c, p, base6, cd, c->r6);
}

static void eightth_order (controls *c, void *p, real h) {
    base8(c, p, h);
}

static void base10 (controls *c, void *p, real cd) {
    suzuki(c, p, base8, cd, c->r8);
}

static void tenth_order (controls *c, void *p, real h) {
    base10(c, p, h);
}

static weights w (int order) {  // composition increases order by one, then symmetry bumps that to the next even order!
    real f = 1.0L / (4.0L - powl(4.0L, 1.0L / (order + 1)));
    return (weights){.fwd = f, .rev = 1.0L - 4.0L * f};
}

static integrator get_integrator (controls *c) {
    integrator composer = NULL;
    switch (c->order) {
        case  2: composer =  second_order; break;
        case  4: composer =  fourth_order; c->r2 = w(2); break;
        case  6: composer =   sixth_order; c->r2 = w(2); c->r4 = w(4); break;
        case  8: composer = eightth_order; c->r2 = w(2); c->r4 = w(4); c->r6 = w(6); break;
        case 10: composer =   tenth_order; c->r2 = w(2); c->r4 = w(4); c->r6 = w(6); c->r8 = w(8); break;
    }
    CHECK(composer);
    return composer;
}

void solve (char **argv, controls *c, void *p, plotter output) {
    int display_precision = (int)strtol(argv[1], NULL, BASE); CHECK(display_precision >= 1 && display_precision <= 32);
    integrator composer = get_integrator(c);
    for (int step = 0; step < c->steps; step++) {
        output(display_precision, p, step * c->step_size);
        composer(c, p, c->step_size);
    }
    output(display_precision, p, c->steps * c->step_size);
}

bool generate (controls *c, void *p) {
    static bool looping = false;
    static integrator composer;
    if (looping) goto resume; else looping = true;
    composer = get_integrator(c);
    for (c->step = 0; c->step < c->steps; c->step++) {
        composer(c, p, c->step_size);
        return true;
        resume: ;
    }
    return looping = false;
}
