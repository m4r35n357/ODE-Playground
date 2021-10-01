/*
 * Wimol-Banlue System
 *
 * Example: ./tsm-wimol-banlue-dbg 9 32 8 0.1 10000 1 0 0 2.0
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a, d0; series tx, s2x; } parameters;

void *get_p (int argc, char **argv, long n) {
    assert(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a);
    mpfr_init(p->d0);
    mpfr_set_ui(p->d0, 0, RND);
    p->tx = t_jet(n); p->s2x = t_jet(n);
    return p;
}

void ode (series x, series y, series z, components *c, void *params, int k) {
    parameters *p = (parameters *)params;
	//  x' = y - x
	mpfr_sub(c->x, y[k], x[k], RND);
	//  y' = - z * tan(x)
	t_tan_sec2(p->tx, p->s2x, x, k, HYP);
    mpfr_neg(c->y, *t_prod(z, p->tx, k), RND);
	//  z' = - A + xy + |y|
	mpfr_add(c->z, *t_prod(x, y, k), *t_abs(y, k), RND);
	mpfr_sub(c->z, c->z, k == 0 ? p->a : p->d0, RND);
}
