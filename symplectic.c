/*
 * Symplectic Integrators and support functions
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symplectic.h"

real error (real e) {
    return - log10l(fabsl(e) >= 1e-36L ? fabsl(e) : 1e-36L);
}

controls *symp_get_c (int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    controls *_ = malloc(sizeof (controls)); CHECK(_);
    _->dp = (int)strtol(argv[1], NULL, BASE);    CHECK(_->dp >= 1);
    _->order = (int)strtol(argv[2], NULL, BASE); CHECK(_->order > 0 && _->order % 2 == 0);
    _->h = strtold(argv[3], NULL);               CHECK(_->h > 0.0L);
    _->steps = (int)strtol(argv[4], NULL, BASE); CHECK(_->steps >= 0 && _->steps <= 1000000);
    _->looping = false;
    return _;
}

static void _symplectic_ (int order, model *p, real c_d) {
    if (order > 2) {
        order -= 2;
        real fwd = 1.0L / (4.0L - powl(4.0L, 1.0L / (order + 1)));
        for (int stage = 0; stage < 5; stage++) {
            _symplectic_(order, p, (stage == 2 ? 1.0L - 4.0L * fwd : fwd) * c_d);
        }
    } else {
        update_q(p, c_d * 0.5L);
        update_p(p, c_d);
        update_q(p, c_d * 0.5L);
    }
}

void solve (controls *c, model *p, plotter output) {
    for (int step = 0; step < c->steps; step++) {
        output(c->dp, p, step * c->h);
        _symplectic_(c->order, p, c->h);
    }
    output(c->dp, p, c->steps * c->h);
}

bool generate (controls *c, model *p) {
    if (c->looping) goto resume; else c->looping = true;
    for (c->step = 0; c->step < c->steps; c->step++) {
        _symplectic_(c->order, p, c->h);
        return true;
        resume: ;
    }
    return c->looping = false;
}
