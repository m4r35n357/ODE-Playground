/*
 * Thomas' cyclically symmetric attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t b; series sx, sy, sz, cx, cy, cz; };

parameters *get_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    tsm_get_p(argv, argc, &p->b);
    p->sx = tsm_jet(n); p->cx = tsm_jet(n);
    p->sy = tsm_jet(n); p->cy = tsm_jet(n);
    p->sz = tsm_jet(n); p->cz = tsm_jet(n);
    return p;
}

void ode (triplet *v_k, series x, series y, series z, parameters *p, int k) {
    //  x' = sin(y) - Bx
    mpfr_fms(v_k->x, p->b, x[k], *t_sin_cos(p->sy, p->cy, y, k, true).a, RND);
    mpfr_neg(v_k->x, v_k->x, RND);
    //  y' = sin(z) - By
    mpfr_fms(v_k->y, p->b, y[k], *t_sin_cos(p->sz, p->cz, z, k, true).a, RND);
    mpfr_neg(v_k->y, v_k->y, RND);
    //  z' = sin(x) - Bz
    mpfr_fms(v_k->z, p->b, z[k], *t_sin_cos(p->sx, p->cx, x, k, true).a, RND);
    mpfr_neg(v_k->z, v_k->z, RND);
}
