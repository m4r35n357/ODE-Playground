/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-dbg  6 8  0.1 30000  1 0 0  .185
 *
 $(yad --columns=2 --title="Thomas Attractor (TSM)" --form --separator=" " --align=right \
    --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="b" \
    -- './tsm-thomas-static!./tsm-thomas-dbg!./tsm-thomas' \
    '6!0..36!3' '8!4..32!1' '.1!0.001..0.1!0.001!3' '30000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '0.185') >/tmp/$USER/data
 *
 ./cns $(yad --columns=2 --title="Thomas CNS (TSM)" --form --separator=" " --align=right \
    --field="Mode":CB --field="Separation" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="b" \
    -- 'step2!nosim' '1.0' './tsm-thomas-static!./tsm-thomas-dbg!./tsm-thomas' \
    '6!0..36!3' '8!4..32!1' '.1!0.001..0.1!0.001!3' '30000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '0.185')
 *
 ./cns-scan $(yad --columns=2 --title="Thomas CNS Scan (TSM)" --form --separator=" " --align=right \
    --field="Maxium Order":NUM --field="Separation" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":RO --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="b" \
    -- '32!2..64!1' '1.0' './tsm-thomas-static!./tsm-thomas-dbg!./tsm-thomas' \
    '6!0..36!3' '_' '.1!0.001..0.1!0.001!3' '30000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '0.185') | tee /tmp/$USER/data
 *
 ./bifurcation-scan $(yad --columns=2 --title="Thomas Bifurcation (TSM)" --form --separator=" " --align=right \
    --field="Lower Value" --field="Upper Value" --field="Skip Transient" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="b:RO" \
    -- '0.1' '0.23' '10' './tsm-thomas-static!./tsm-thomas-dbg!./tsm-thomas' \
    '6!0..36!3' '4!4..32!1' '.1!0.001..0.1!0.001!3' '30000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '$p')
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real b; series sx, sy, sz, cx, cy, cz; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 9);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->b);
    p->sx = t_jet(n); p->cx = t_jet(n);
    p->sy = t_jet(n); p->cy = t_jet(n);
    p->sz = t_jet(n); p->cz = t_jet(n);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = t_sin_cos(p->sy, p->cy, y, k, TRIG).a - p->b * x[k],
        .y = t_sin_cos(p->sz, p->cz, z, k, TRIG).a - p->b * y[k],
        .z = t_sin_cos(p->sx, p->cx, x, k, TRIG).a - p->b * z[k]
    };
}
