/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-dbg  6 8  .01 10000  1 0 0  1.4 4
 *
 $(yad --columns=2 --title="Halvorsen Attractor (TSM)" --form --separator=" " --align=right \
    --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="a" --field="b" \
    -- './tsm-halvorsen-static!./tsm-halvorsen-dbg!./tsm-halvorsen!./tsm-halvorsen.py' \
    '6!3..64!3' '8!4..128!1' '.1!0.001..0.01!0.001!3' '10000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '1.4' '4') >/tmp/$USER/data
 *
 ./cns $(yad --columns=2 --title="Halvorsen CNS (TSM)" --form --separator=" " --align=right \
    --field="Mode":CB --field="Separation" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="a" --field="b" \
    -- 'step2!nosim' '1.0' './tsm-halvorsen-static!./tsm-halvorsen-dbg!./tsm-halvorsen!./tsm-halvorsen.py' \
    '6!3..64!3' '8!4..128!1' '.1!0.001..0.01!0.001!3' '10000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '1.4' '4')
 *
 ./cns-scan $(yad --columns=2 --title="Halvorsen CNS Scan (TSM)" --form --separator=" " --align=right \
    --field="Maxium Order":NUM --field="Separation" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":RO --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="a" --field="b" \
    -- '32!2..64!1' '1.0' './tsm-halvorsen-static!./tsm-halvorsen-dbg!./tsm-halvorsen!./tsm-halvorsen.py' \
    '6!3..64!3' '_' '.1!0.001..0.01!0.001!3' '10000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '1.4' '4') | tee /tmp/$USER/data
 *
 ./bifurcation-scan $(yad --title="Halvorsen Bifurcation (TSM)" --form --separator=" " --align=right \
    --field="Lower Value" --field="Upper Value" --field="Skip Transient" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="a:RO" --field="b" \
    -- '0.1' '0.23' '10' './tsm-halvorsen-static!./tsm-halvorsen-dbg!./tsm-halvorsen!./tsm-halvorsen.py' \
    '6!3..64!3' '4!4..128!1' '.1!0.001..0.01!0.001!3' '10000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '$p' '4')
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real a, b; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 10);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = - p->a * x[k] - p->b * (y[k] + z[k]) - t_sqr(y, k),
        .y = - p->a * y[k] - p->b * (z[k] + x[k]) - t_sqr(z, k),
        .z = - p->a * z[k] - p->b * (x[k] + y[k]) - t_sqr(x, k)
    };
}
