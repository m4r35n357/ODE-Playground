/*
 * Inverted smooth unimodal chaos - http://www.atomosyd.net/spip.php?article218
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, b, c; series x2py2; };

model *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 11);
    model *_ = malloc(sizeof (model)); CHECK(_);
    _->x2py2 = tsm_jet(n);
    tsm_get_p(argv, argc, &_->a, &_->b, &_->c);
    return _;
}

triplet ode (series x, series y, series z, const model *_, int k) {
    _->x2py2[k] = t_sqr(x, k) + t_sqr(y, k);
    return (triplet) {
        .x = z[k] - y[k],
        .y = x[k] - _->a * y[k],
        .z = t_const(_->b, k) + _->c * z[k] - t_mul(_->x2py2, z, k)
    };
}
