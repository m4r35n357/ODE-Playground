/*
 * Rossler System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t a, b, c; series _B; };

parameters *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 12);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c);
    p->_B = t_const(n, p->b);
    return p;
}

void ode (triplet *v_k, series x, series y, series z, parameters *p, int k) {
    //  x' = - y - z
    mpfr_add(v_k->x, y[k], z[k], RND);
    mpfr_neg(v_k->x, v_k->x, RND);
    //  y' = x + Ay
    mpfr_fma(v_k->y, p->a, y[k], x[k], RND);
    //  z' = B + z(x - C)
    mpfr_fms(v_k->z, p->c, z[k], *t_mul(z, x, k), RND);
    mpfr_sub(v_k->z, p->_B[k], v_k->z, RND);
}
