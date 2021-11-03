/*
 * Burke & Shaw System http://www.atomosyd.net/spip.php?article33
 *
 * Example: ./tsm-burke-shaw-dbg 9 10 .01 10000 1 1 1 10 4.272
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real s, v; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 10);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->s, &p->v);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = - p->s * (x[k] + y[k]),
        .y = - (p->s * t_prod(x, z, k) + y[k]),
        .z = p->s * t_prod(x, y, k) + p->v
    };
}
