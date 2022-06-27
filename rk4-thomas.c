/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./rk4-thomas-dbg  15 1  0.1 30000  1 0 0  .185
 *
 $(yad --columns=2 --title="Thomas Attractor (RK4)" --form --separator=" " --align=right \
    --field="Model:CB" \
    --field="Display Places":NUM --field="Plot Interval":NUM --field="Step Size" --field="Steps" \
    --field="x0" --field="y0" --field="z0" \
    --field="b" \
    -- './rk4-thomas-static!./rk4-thomas-dbg!./rk4-thomas' \
    '6!0..36!3' '1!1..1000!1' '.1' '30000' \
    '1.0' '0.0' '0.0' \
    '0.185') >/tmp/$USER/data
 *
 ./cns $(yad --columns=2 --title="Thomas CNS (RK4)" --form --separator=" " --align=right \
    --field="Mode":CB --field="Separation" --field="Model:CB" \
    --field="Display Places":NUM --field="Plot Interval":NUM --field="Step Size" --field="Steps" \
    --field="x0" --field="y0" --field="z0" \
    --field="b" \
    -- 'step2!nosim' '1.0' './rk4-thomas-static!./rk4-thomas-dbg!./rk4-thomas' \
    '6!0..36!3' '8!4..128!1' '.1' '30000' \
    '1.0' '0.0' '0.0' \
    '0.185')
 *
 ./bifurcation-scan $(yad --columns=2 --title="Thomas Bifurcation (RK4)" --form --separator=" " --align=right \
    --field="Lower Value" --field="Upper Value" --field="Skip Transient" --field="Model:CB" \
    --field="Display Places":NUM --field="Plot Interval":NUM --field="Step Size" --field="Steps" \
    --field="x0" --field="y0" --field="z0" \
    --field="b:RO" \
    -- '0.1' '0.23' '10' './rk4-thomas-static!./rk4-thomas-dbg!./rk4-thomas' \
    '6!0..36!3' '4!4..128!1' '.1' '30000' \
    '1.0' '0.0' '0.0' \
    '$p')
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "runge-kutta-ode.h"

typedef struct { real b; } parameters;

void *get_p (int argc, char **argv) {
    assert(argc == 9);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->b);
    return p;
}

components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = sinl(y) - p->b * x,
        .y = sinl(z) - p->b * y,
        .z = sinl(x) - p->b * z
    };
}
