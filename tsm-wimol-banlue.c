/*
 * Wimol-Banlue System
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t a; series tx, s2x, _A; };

parameters *get_p (int argc, char **argv, int n) {
    CHECK(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a);
    p->tx = t_jet(n);
    p->s2x = t_jet(n);
    p->_A = t_const(n, p->a);
    return p;
}

void ode (components *vk, series x, series y, series z, parameters *p, int k) {
    //  x' = y - x
    mpfr_sub(vk->x, y[k], x[k], RND);
    //  y' = - z * tan(x)
    t_tan_sec2(p->tx, p->s2x, x, k, false);
    mpfr_neg(vk->y, *t_mul(z, p->tx, k), RND);
    //  z' = - A + xy + |y|
    mpfr_add(vk->z, *t_mul(x, y, k), *t_abs(y, k), RND);
    mpfr_sub(vk->z, vk->z, p->_A[k], RND);
}
