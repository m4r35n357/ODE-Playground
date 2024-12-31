/*
 * Inverted smooth unimodal chaos http://www.atomosyd.net/spip.php?article218
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "taylor-ode.h"

struct Parameters { real a, b, c; series x2py2, _B; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 12);
    model *_ = malloc(sizeof (model)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b, &_->c);
    _->x2py2 = tsm_jet(n);
    _->_B = tsm_jet(n); mpfr_set(_->_B[0], _->b, RND);
    return _;
}

void ode (triplet *v, series x, series y, series z, model *_, int k) {
    //  x' = z - y
    mpfr_sub(v->x, z[k], y[k], RND);
    //  y' = x - Ay
    mpfr_fms(v->y, _->a, y[k], x[k], RND);
    mpfr_neg(v->y, v->y, RND);
    //  z' = B + Cz - (x^2 + y^2)z
    mpfr_add(_->x2py2[k], *t_sqr(x, k), *t_mul(y, y, k), RND);
    mpfr_fms(v->z, _->c, z[k], *t_mul(_->x2py2, z, k), RND);
    mpfr_add(v->z, v->z, _->_B[k], RND);
}
