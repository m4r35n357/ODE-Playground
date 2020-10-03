/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./rk4-rf-dbg 15 NA 1 .01 50000 .05 -.05 .3 .28713 .1
 *          ./rk4-rf-dbg 15 NA 1 .01 50000 .05 -.05 .3 .105 .1 | ./plotPi3d.py
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct {
    real alpha;
    real gamma;
} parameters;

static void *get_p (int argc, char **argv) {
    parameters *p = malloc(sizeof (parameters));
    t_args(argv, argc, &p->alpha, &p->gamma);
    return p;
}

static components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = y * (z - 1.0 + x * x) + p->gamma * x,
        .y = x * (3.0 * z + 1.0 - x * x) + p->gamma * y,
        .z = - 2.0 * z * (p->alpha + x * y)
    };
}

int main (int argc, char **argv) {
    assert(argc == 11);
    rk4(argc, argv, ode, get_p);
    return 0;
}
