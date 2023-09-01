/*
 * Yu-Wang System
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t a, b, c, d; series xy, e_xy; };

parameters *get_p (int argc, char **argv, int n) {
    CHECK(argc == 13);
    parameters *p = malloc(sizeof (parameters));
    tsm_get_p(argv, argc, &p->a, &p->b, &p->c, &p->d);
    p->xy = tsm_var(n); p->e_xy = tsm_var(n);
    return p;
}

void ode (triplet *v_k, series x, series y, series z, parameters *p, int k) {
    //  x' = A(y - x)
    mpfr_fmms(v_k->x, p->a, y[k], p->a, x[k], RND);
    //  y' = Bx - cxz
    mpfr_fmms(v_k->y, p->b, x[k], p->c, *t_mul(x, z, k), RND);
    //  z' = e^(xy) - Dz
    mpfr_swap(p->xy[k], *t_mul(x, y, k));
    mpfr_fms(v_k->z, p->d, z[k], *t_exp(p->e_xy, p->xy, k), RND);
    mpfr_neg(v_k->z, v_k->z, RND);
}
