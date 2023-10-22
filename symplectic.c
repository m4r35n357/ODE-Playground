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
    _->dp = (int)strtol(argv[1], NULL, BASE);    //CHECK(_->dp >= 1 && _->dp <= 32);
    _->order = (int)strtol(argv[2], NULL, BASE); CHECK(_->order > 0 && _->order % 2 == 0);
    _->h = strtold(argv[3], NULL);               CHECK(_->h > 0.0L);
    _->steps = (int)strtol(argv[4], NULL, BASE); CHECK(_->steps >= 0 && _->steps <= 1000000);
    _->looping = false;
    return _;
}

real error (real e) {
    return - log10l(fabsl(e) >= 1e-36L ? fabsl(e) : 1e-36L);
}

static void suzuki (controls *c, model *p, block base, real cd, pair weight) {
    for (int stage = 0; stage < 5; stage++) {
        base(c, p, (stage == 2 ? weight.b : weight.a) * cd);
    }
}

static void _base2_ (controls *cont, model *p, real cd) { (void)cont;
    update_q(p, cd * 0.5L);
    update_p(p, cd);
    update_q(p, cd * 0.5L);
}

static void _base4_ (controls *c, model *p, real cd) {
    suzuki(c, p, _base2_, cd, c->w2);
}

static void _base6_ (controls *c, model *p, real cd) {
    suzuki(c, p, _base4_, cd, c->w4);
}

static void _base8_ (controls *c, model *p, real cd) {
    suzuki(c, p, _base6_, cd, c->w6);
}

static void _base10_ (controls *c, model *p, real cd) {
    suzuki(c, p, _base8_, cd, c->w8);
}

static void second_order (controls *c, model *p) {
    _base2_(c, p, c->h);
}

static void fourth_order (controls *c, model *p) {
    _base4_(c, p, c->h);
}

static void sixth_order (controls *c, model *p) {
    _base6_(c, p, c->h);
}

static void eightth_order (controls *c, model *p) {
    _base8_(c, p, c->h);
}

static void tenth_order (controls *c, model *p) {
    _base10_(c, p, c->h);
}

static pair w (int order) {  // composition increases order by one, then symmetry bumps that to the next even order!
    real _ = 1.0L / (4.0L - powl(4.0L, 1.0L / (order + 1)));
    return (pair){.a = _, .b = 1.0L - 4.0L * _};
}

static integrator get_integrator (controls *c) {
    integrator _ = NULL;
    switch (c->order) {
        case  2: _ =  second_order; break;
        case  4: _ =  fourth_order; c->w2 = w(2); break;
        case  6: _ =   sixth_order; c->w2 = w(2); c->w4 = w(4); break;
        case  8: _ = eightth_order; c->w2 = w(2); c->w4 = w(4); c->w6 = w(6); break;
        case 10: _ =   tenth_order; c->w2 = w(2); c->w4 = w(4); c->w6 = w(6); c->w8 = w(8); break;
        default: break;
    }; CHECK(_);
    return _;
}

void solve (controls *c, model *p, plotter output) {
    integrator symplectic = get_integrator(c);
    for (int step = 0; step < c->steps; step++) {
        output(c->dp, p, step * c->h);
        symplectic(c, p);
    }
    output(c->dp, p, c->steps * c->h);
}

bool generate (controls *c, model *p) {
    integrator symplectic = get_integrator(c);
    if (c->looping) goto resume; else c->looping = true;
    for (c->step = 0; c->step < c->steps; c->step++) {
        symplectic(c, p);
        return true;
        resume: ;
    }
    return c->looping = false;
}
