/*
 * Burke & Shaw System http://www.atomosyd.net/spip.php?article33
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t s, v; };

parameters *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 11);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->s, &p->v);
    return p;
}

void ode (components *v_k, series x, series y, series z, parameters *p, int k) {
    //  x' = - S(x + y)
    mpfr_fmma(v_k->x, p->s, x[k], p->s, y[k], RND);
    mpfr_neg(v_k->x, v_k->x, RND);
    //  y' = - (Sxz + y)
    mpfr_fma(v_k->y, p->s, *t_mul(x, z, k), y[k], RND);
    mpfr_neg(v_k->y, v_k->y, RND);
    //  z' = Sxy + V
    mpfr_fma(v_k->z, p->s, *t_mul(x, y, k), p->v, RND);
}
