/*
 * Burke & Shaw System http://www.atomosyd.net/spip.php?article33
 *
 * Example: ./tsm-burke-shaw-dbg 9 32 10 .01 10000 1 1 1 10 4.272
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t s, v; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 11);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->s, &p->v);
    return p;
}

void ode (components *vk, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = - S(x + y)
    mpfr_fmma(vk->x, p->s, x[k], p->s, y[k], RND);
    mpfr_neg(vk->x, vk->x, RND);
    //  y' = - (Sxz + y)
    mpfr_fma(vk->y, p->s, *t_mul(x, z, k), y[k], RND);
    mpfr_neg(vk->y, vk->y, RND);
    //  z' = Sxy + V
    mpfr_fma(vk->z, p->s, *t_mul(x, y, k), p->v, RND);
}
