/*
 * Rucklidge Attractor
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { mpfr_t alpha, kappa; };

parameters *get_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 11);
    parameters *p = malloc(sizeof (parameters));
    tsm_get_p(argv, argc, &p->alpha, &p->kappa);
    return p;
}

void ode (triplet *v_k, series x, series y, series z, parameters *p, int k) {
    //  x' = ay - kx - yz
    mpfr_fmms(v_k->x, p->alpha, y[k], p->kappa, x[k], RND);
    mpfr_sub(v_k->x, v_k->x, *t_mul(y, z, k), RND);
    //  y' = x
    mpfr_set(v_k->y, x[k], RND);
    //  z' = y^2 - z
    mpfr_sub(v_k->z, *t_sqr(y, k), z[k], RND);
}
