/*
 * Genesio-Tesi System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
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

void ode (triplet *v_k, series x, series y, series z, parameters *p, int k) {
    //  x' = y
    mpfr_set(v_k->x, y[k], RND);
    //  y' = z
    mpfr_set(v_k->y, z[k], RND);
    //  z' = - Az - By - x(1 + x)
    mpfr_fmma(v_k->z, p->a, z[k], p->b, y[k], RND);
    mpfr_add(v_k->z, v_k->z, x[k], RND);
    mpfr_add(v_k->z, v_k->z, *t_sqr(x, k), RND);
    mpfr_neg(v_k->z, v_k->z, RND);
}
