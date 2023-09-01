/*
 * Halvorsen Cyclic Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t a, D4; };

parameters *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    tsm_get_p(argv, argc, &p->a);
    mpfr_init_set_si(p->D4, 4, RND);
    return p;
}

void ode (triplet *v_k, series x, series y, series z, parameters *p, int k) {
    //  x' = - Ax - 4y - 4z - y^2
    mpfr_fmma(v_k->x, p->D4, y[k], p->D4, z[k], RND);
    mpfr_fma(v_k->x, p->a, x[k], v_k->x, RND);
    mpfr_add(v_k->x, *t_sqr(y, k), v_k->x, RND);
    mpfr_neg(v_k->x, v_k->x, RND);
    //  y' = - Ay - 4z - 4x - z^2
    mpfr_fmma(v_k->y, p->D4, z[k], p->D4, x[k], RND);
    mpfr_fma(v_k->y, p->a, y[k], v_k->y, RND);
    mpfr_add(v_k->y, *t_sqr(z, k), v_k->y, RND);
    mpfr_neg(v_k->y, v_k->y, RND);
    //  z' = - Az - 4x - 4y - x^2
    mpfr_fmma(v_k->z, p->D4, x[k], p->D4, y[k], RND);
    mpfr_fma(v_k->z, p->a, z[k], v_k->z, RND);
    mpfr_add(v_k->z, *t_sqr(x, k), v_k->z, RND);
    mpfr_neg(v_k->z, v_k->z, RND);
}
