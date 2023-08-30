/*
 * Inverted smooth unimodal chaos http://www.atomosyd.net/spip.php?article218
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a, b, c; series x2py2, _B; };

parameters *tsm_init_p (int argc, char **argv, int n) {
    CHECK(argc == 11);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->a, &_->b, &_->c);
    _->x2py2 = tsm_var(n);
    _->_B = tsm_const(n, _->b);
    return _;
}

triplet ode (series x, series y, series z, parameters *_, int k) {
    _->x2py2[k] = t_sqr(x, k) + t_sqr(y, k);
    return (triplet) {
        .x = z[k] - y[k],
        .y = x[k] - _->a * y[k],
        .z = _->_B[k] + _->c * z[k] - t_mul(_->x2py2, z, k)
    };
}
