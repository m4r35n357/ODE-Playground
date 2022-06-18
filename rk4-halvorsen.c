/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./rk4-halvorsen-dbg  6 1  .01 10000  1 0 0  1.4 4
 *
 $(yad --columns=2 --title="Halvorsen Attractor (RK4)" --form --separator=" " --align=right \
    --field="Model:CB" \
    --field="Display Places":NUM --field="Plot Interval":NUM --field="Step Size" --field="Steps" \
    --field="x0" --field="y0" --field="z0" \
    --field="a" --field="b" \
    -- './rk4-halvorsen-static!./rk4-halvorsen-dbg!./rk4-halvorsen!./rk4-halvorsen.py' \
    '6!3..64!3' '1!1..1000!1' '.01' '10000' \
    '1.0' '0.0' '0.0' \
    '1.4' '4') >/tmp/$USER/data
 *
 ./cns $(yad --columns=2 --title="Halvorsen CNS (RK4)" --form --separator=" " --align=right \
    --field="Mode":CB --field="Separation" --field="Model:CB" \
    --field="Display Places":NUM --field="Plot Interval":NUM --field="Step Size" --field="Steps" \
    --field="x0" --field="y0" --field="z0" \
    --field="a" --field="b" \
    -- 'step2!nosim' '1.0' './rk4-halvorsen-static!./rk4-halvorsen-dbg!./rk4-halvorsen!./rk4-halvorsen.py' \
    '6!3..64!3' '1!1..1000!1' '.01' '10000' \
    '1.0' '0.0' '0.0' \
    '1.4' '4')
 *
 ./bifurcation-scan $(yad --columns=2 --title="Halvorsen Bifurcation (RK4)" --form --separator=" " --align=right \
    --field="Lower Value" --field="Upper Value" --field="Skip Transient" --field="Model:CB" \
    --field="Display Places":NUM --field="Plot Interval":NUM --field="Step Size" --field="Steps" \
    --field="x0" --field="y0" --field="z0" \
    --field="a:RO" --field="b" \
    -- '0.1' '0.23' '10' './rk4-halvorsen-static!./rk4-halvorsen-dbg!./rk4-halvorsen' \
    '6!3..64!3' '1!1..1000!1' '.01' '10000' \
    '1.0' '0.0' '0.0' \
    '$p' '4')
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "runge-kutta-ode.h"

typedef struct { real a, b; } parameters;

void *get_p (int argc, char **argv) {
    assert(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    return p;
}

components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = - p->a * x - p->b * (y + z) - y * y,
        .y = - p->a * y - p->b * (z + x) - z * z,
        .z = - p->a * z - p->b * (x + y) - x * x
    };
}
