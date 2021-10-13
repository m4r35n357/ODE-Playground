/*
 * Logistic Equation
 *
 * Example:  ./tsm-logistic-dbg 9 32 10 0.005 10000 .001 0 0 .5
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 10);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a);
    return p;
}

void ode (components *v, series x, series y, series z, void *params, int k) {
    (void)y; (void)z;
    parameters *p = (parameters *)params;
    //  x' = Ax(1 - x)
    mpfr_fmms(v->x, p->a, x[k], p->a, *t_sqr(x, k), RND);
    // not used
    mpfr_set_zero(v->y, 1);
    // not used
    mpfr_set_zero(v->z, 1);
}
