/*
 * Bouali Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t a, b, g, m; series gx2; };

parameters *get_p (int argc, char **argv, int n) {
    CHECK(argc == 13);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->g, &p->m);
    p->gx2 = t_jet(n);
    return p;
}

void ode (triplet *v_k, series x, series y, series z, parameters *p, int k) {
    //  x' = Ax(1 - y) - Bz
    mpfr_fmms(v_k->x, p->a, x[k], p->a, *t_mul(x, y, k), RND);
    mpfr_fms(v_k->x, p->b, z[k], v_k->x, RND);
    mpfr_neg(v_k->x, v_k->x, RND);
    //  y' = - Gy(1 - x^2)
    mpfr_mul(p->gx2[k], p->g, *t_sqr(x, k), RND);
    mpfr_fms(v_k->y, p->g, y[k], *t_mul(y, p->gx2, k), RND);
    mpfr_neg(v_k->y, v_k->y, RND);
    //  z' = Mx
    mpfr_mul(v_k->z, p->m, x[k], RND);
}
