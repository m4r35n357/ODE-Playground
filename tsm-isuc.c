/*
 * Inverted smooth unimodal chaos http://www.atomosyd.net/spip.php?article218
 *
 * Example: ./tsm-isuc-dbg 9 10 .01 10000 .05 -.05 .3 .1 .1 4.2
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct Parameters { real a, b, c; series x2py2; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 11);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b, &p->c);
    p->x2py2 = t_jet(n);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->x2py2[k] = t_sqr(x, k) + t_sqr(y, k);
    return (components) {
        .x = z[k] - y[k],
        .y = x[k] - p->a * y[k],
        .z = t_const(p->b, k) + p->c * z[k] - t_mul(p->x2py2, z, k)
    };
}
