/*
 * Nose-Hoover System (Sprott A) - http://www.atomosyd.net/spip.php?article193
 *
 * Example: ./tsm-nose-hoover-std 15 10 0.01 10000 1 0 0 6.0
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "taylor-ode.h"

struct Parameters { real a; };

parameters *tsm_init_p (int argc, char **argv, int n) { (void)n;
    CHECK(argc == 9);
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    tsm_get_p(argv, argc, &_->a);
    return _;
}

triplet ode (series x, series y, series z, parameters *p, int k) {
    return (triplet) {
        .x = y[k],
        .y = t_mul(y, z, k) - x[k],
        .z = p->a - t_sqr(y, k)
    };
}
