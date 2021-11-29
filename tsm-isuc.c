/*
 * Inverted smooth unimodal chaos http://www.atomosyd.net/spip.php?article218
 *
 * Example: ./tsm-isuc-dbg 9 32 10 .01 10000 .05 -.05 .3 .1 .1 4.2
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a, b, c; series x2py2; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 12);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c);
    p->x2py2 = t_jet(n);
    return p;
}

void ode (components *v, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = z - y
    mpfr_sub(v->x, z[k], y[k], RND);
    //  y' = x - Ay
    mpfr_fms(v->y, p->a, y[k], x[k], RND);
    mpfr_neg(v->y, v->y, RND);
    //  z' = B + Cz - (x^2 + y^2)z
    mpfr_add(p->x2py2[k], *t_sqr(x, k), *t_mul(y, y, k), RND);
    mpfr_fms(v->z, p->c, z[k], *t_mul(p->x2py2, z, k), RND);
    mpfr_add(v->z, v->z, *t_const(p->b, k), RND);
}
