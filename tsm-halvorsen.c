/*
 * Halvorsen Cyclic Attractor
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t a, b; };

parameters *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 11);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    return p;
}

void ode (components *v_k, series x, series y, series z, parameters *p, int k) {
    //  x' = - Ax - By - Bz - y^2
    mpfr_fmma(v_k->x, p->b, y[k], p->b, z[k], RND);
    mpfr_fma(v_k->x, p->a, x[k], v_k->x, RND);
    mpfr_add(v_k->x, *t_sqr(y, k), v_k->x, RND);
    mpfr_neg(v_k->x, v_k->x, RND);
    //  y' = - Ay - Bz - Bx - z^2
    mpfr_fmma(v_k->y, p->b, z[k], p->b, x[k], RND);
    mpfr_fma(v_k->y, p->a, y[k], v_k->y, RND);
    mpfr_add(v_k->y, *t_sqr(z, k), v_k->y, RND);
    mpfr_neg(v_k->y, v_k->y, RND);
    //  z' = - Az - Bx - By - x^2
    mpfr_fmma(v_k->z, p->b, x[k], p->b, y[k], RND);
    mpfr_fma(v_k->z, p->a, z[k], v_k->z, RND);
    mpfr_add(v_k->z, *t_sqr(x, k), v_k->z, RND);
    mpfr_neg(v_k->z, v_k->z, RND);
}
