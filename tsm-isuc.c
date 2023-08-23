/*
 * Inverted smooth unimodal chaos http://www.atomosyd.net/spip.php?article218
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t a, b, c; series x2py2, _B; };

parameters *get_p (int argc, char **argv, int n) {
    CHECK(argc == 12);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c);
    p->x2py2 = t_jet(n);
    p->_B = t_const(n, p->b);
    return p;
}

void ode (triplet *v_k, series x, series y, series z, parameters *p, int k) {
    //  x' = z - y
    mpfr_sub(v_k->x, z[k], y[k], RND);
    //  y' = x - Ay
    mpfr_fms(v_k->y, p->a, y[k], x[k], RND);
    mpfr_neg(v_k->y, v_k->y, RND);
    //  z' = B + Cz - (x^2 + y^2)z
    mpfr_add(p->x2py2[k], *t_sqr(x, k), *t_mul(y, y, k), RND);
    mpfr_fms(v_k->z, p->c, z[k], *t_mul(p->x2py2, z, k), RND);
    mpfr_add(v_k->z, v_k->z, p->_B[k], RND);
}
