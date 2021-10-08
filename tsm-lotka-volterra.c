/*
 * Lotka-Volterra (Predator-Prey) System
 *
 * Example: ./tsm-lotka-volterra-dbg 9 32 10 .01 2000 10 10 0 1 .5 .05 .02
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a, b, c, d, xy; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 13);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c, &p->d);
    mpfr_init(p->xy);
    return p;
}

void ode (series x, series y, series z, components *c, void *params, int k) {
    (void)z;
    parameters *p = (parameters *)params;
    mpfr_swap(p->xy, *t_prod(x, y, k));
    //  x' = Ax - Cxy
    mpfr_fmms(c->x, p->a, x[k], p->c, p->xy, RND);
    //  y' = Dxy - By
    mpfr_fmms(c->y, p->d, p->xy, p->b, y[k], RND);
    // not used
    mpfr_set_zero(c->z, 1);
}
