/*
 * Bouali Attractor
 *
 * Example: ./tsm-bouali-dbg 9 32 10 0.01 50000 1 1 0 3 2.2 1 .01
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a, b, g, m; series gx2; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 13);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->g, &p->m);
    p->gx2 = t_jet(n);
    return p;
}

void ode (components *v, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = Ax(1 - y) - Bz
    mpfr_fmms(v->x, p->a, x[k], p->a, *t_mul(x, y, k), RND);
    mpfr_fms(v->x, p->b, z[k], v->x, RND);
    mpfr_neg(v->x, v->x, RND);
    //  y' = - Gy(1 - x^2)
    mpfr_mul(p->gx2[k], p->g, *t_sqr(x, k), RND);
    mpfr_fms(v->y, p->g, y[k], *t_mul(y, p->gx2, k), RND);
    mpfr_neg(v->y, v->y, RND);
    //  z' = Mx
    mpfr_mul(v->z, p->m, x[k], RND);
}
