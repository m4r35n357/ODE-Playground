/*
 * Logistic Equation
 *
 * Example:  ./tsm-logistic-std 9 10 0.005 10000 .001 0 0 .5
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct Parameters { real a; } parameters;

void *get_p (int argc, char **argv, int n) { (void)n;
    assert(argc == 9);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a);
    return p;
}

components ode (series x, series y, series z, void *params, int k) { (void)y; (void)z;
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->a * (x[k] - t_sqr(x, k))
    };
}
