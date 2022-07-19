/*
 * Rabinovich–Fabrikant System
 *
 * Example: ./tsm-rf-dbg 15 10 .01 50000 .05 -.05 .3 .2873 .1
 * Example: ./tsm-rf-dbg 15 12 .01 100000 .05 -.05 .3 .116364 .1
 * Example: ./tsm-rf-dbg 15 16 .01 50000 .05 -.05 .3 .105 .1
 *
 $(yad --columns=2 --title="Rabinovich–Fabrikant Attractor (TSM)" --form --separator=" " --align=right \
    --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="alpha" --field="gamma" \
    -- './tsm-rf-static!./tsm-rf-dbg!./tsm-rf!./tsm-rf-gl' \
    '6!0..36!3' '8!4..32!1' '.01!0.001..0.1!0.001!3' '50000!1..1000000!1000' \
    '0.05' '-0.05' '0.3' \
    '.2873' '.1') >/tmp/$USER/data
 *
 ./cns $(yad --columns=2 --title="Rabinovich–Fabrikant CNS (TSM)" --form --separator=" " --align=right \
    --field="Mode":CB --field="Separation" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="alpha" --field="gamma" \
    -- 'step2!nosim' '1.0' './tsm-rf-static!./tsm-rf-dbg!./tsm-rf!./tsm-rf-gl' \
    '6!0..36!3' '8!4..32!1' '0.1!0.001..0.1!0.001!3' '30000!1..1000000!1000' \
    '0.05' '-0.05' '0.3' \
    '.2873' '.1')
 *
 ./cns-scan $(yad --columns=2 --title="Rabinovich–Fabrikant CNS Scan (TSM)" --form --separator=" " --align=right \
    --field="Maxium Order":NUM --field="Separation" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":RO --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="alpha" --field="gamma" \
    -- '32!2..64!1' '1.0' './tsm-rf-static!./tsm-rf-dbg!./tsm-rf' \
    '6!0..36!3' '_' '.01!0.001..0.1!0.001!3' '30000!1..1000000!1000' \
    '0.05' '-0.05' '0.3' \
    '.2873' '.1') | tee /tmp/$USER/data
 *
 ./bifurcation-scan $(yad --columns=2 --title="Rabinovich–Fabrikant Bifurcation (TSM)" --form --separator=" " --align=right \
    --field="Lower Value" --field="Upper Value" --field="Skip Transient" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="alpha":RO --field="gamma" \
    -- '0.1' '0.23' '10' './tsm-rf-static!./tsm-rf-dbg!./tsm-rf' \
    '6!0..36!3' '4!4..32!1' '.01!0.001..0.1!0.001!3' '30000!1..1000000!1000' \
    '0.05' '-0.05' '0.3' \
    '$p' '.1')
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct Parameters { real alpha, gamma; series a, b, c; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->alpha, &p->gamma);
    p->a = t_jet(n);
    p->b = t_jet(n);
    p->c = t_jet(n);
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    p->a[k] = z[k] + t_sqr(x, k) - t_const(1.0L, k);
    p->b[k] = 4.0L * z[k] - p->a[k];
    p->c[k] = t_const(p->alpha, k) + t_mul(x, y, k);
    return (components) {
        .x = t_mul(y, p->a, k) + p->gamma * x[k],
        .y = t_mul(x, p->b, k) + p->gamma * y[k],
        .z = - 2.0L * t_mul(z, p->c, k)
    };
}
