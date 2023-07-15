/*
 * Symplectic Integrators and support functions
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symplectic.h"

controls *symp_get_c (int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    controls *_ = malloc(sizeof (controls)); CHECK(_);
    _->order = (int)strtol(argv[2], NULL, BASE); CHECK(_->order > 0 && _->order % 2 == 0);
    _->step_size = strtold(argv[3], NULL);       CHECK(_->step_size > 0.0L);
    _->steps = (int)strtol(argv[4], NULL, BASE); CHECK(_->steps >= 0 && _->steps <= 1000000);
    return _;
}

real error (real e) {
    return - log10l(fabsl(e) >= 1e-36L ? fabsl(e) : 1e-36L);
}

static void second_order (controls *cont, parameters *p, real cd) { (void)cont;
    update_q(p, cd * 0.5L);
    update_p(p, cd);
    update_q(p, cd * 0.5L);
}

static void suzuki (controls *c, parameters *p, integrator base, real cd, pair weights) {
    for (int step = 0; step < 5; step++) base(c, p, cd * (step == 2 ? weights.b : weights.a));
}

static void _base4 (controls *c, parameters *p, real cd) {
    suzuki(c, p, second_order, cd, c->r2);
}

static void fourth_order (controls *c, parameters *p, real h) {
    _base4(c, p, h);
}

static void _base6 (controls *c, parameters *p, real cd) {
    suzuki(c, p, _base4, cd, c->r4);
}

static void sixth_order (controls *c, parameters *p, real h) {
    _base6(c, p, h);
}

static void _base8 (controls *c, parameters *p, real cd) {
    suzuki(c, p, _base6, cd, c->r6);
}

static void eightth_order (controls *c, parameters *p, real h) {
    _base8(c, p, h);
}

static void _base10 (controls *c, parameters *p, real cd) {
    suzuki(c, p, _base8, cd, c->r8);
}

static void tenth_order (controls *c, parameters *p, real h) {
    _base10(c, p, h);
}

static pair w (int order) {  // composition increases order by one, then symmetry bumps that to the next even order!
    real _ = 1.0L / (4.0L - powl(4.0L, 1.0L / (order + 1)));
    return (pair){.a = _, .b = 1.0L - 4.0L * _};
}

static integrator get_integrator (controls *c) {
    integrator _ = NULL;
    switch (c->order) {
        case  2: _ =  second_order; break;
        case  4: _ =  fourth_order; c->r2 = w(2); break;
        case  6: _ =   sixth_order; c->r2 = w(2); c->r4 = w(4); break;
        case  8: _ = eightth_order; c->r2 = w(2); c->r4 = w(4); c->r6 = w(6); break;
        case 10: _ =   tenth_order; c->r2 = w(2); c->r4 = w(4); c->r6 = w(6); c->r8 = w(8); break;
        default: break;
    }; CHECK(_);
    return _;
}

void solve (char **argv, controls *c, parameters *p, plotter output) {
    int _ = (int)strtol(argv[1], NULL, BASE); CHECK(_ >= 1 && _ <= 32); // display precision
    integrator evolve = get_integrator(c);
    for (int step = 0; step < c->steps; step++) {
        output(_, p, step * c->step_size);
        evolve(c, p, c->step_size);
    }
    output(_, p, c->steps * c->step_size);
}

bool generate (controls *c, parameters *p) {
    static bool looping = false;
    static integrator evolve;
    if (looping) goto resume; else looping = true;
    evolve = get_integrator(c);
    for (c->step = 0; c->step < c->steps; c->step++) {
        evolve(c, p, c->step_size);
        return true;
        resume: ;
    }
    return looping = false;
}
