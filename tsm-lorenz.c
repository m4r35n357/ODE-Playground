/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-dbg  6 8 .01 10000  -15.8 -17.48 35.64  10 28 8 3
 *
 $(yad --columns=2 --title="Lorenz Attractor (TSM)" --form --separator=" " --align=right \
    --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="sigma" --field="rho" --field="beta (numerator)" --field="beta (denominator)" \
    -- './tsm-lorenz-static!./tsm-lorenz-dbg!./tsm-lorenz!./tsm-lorenz.py' \
    '6!3..64!3' '8!4..128!1' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' \
    '-15.8' '-18.48' '35.64' \
    '10' '28' '8' '3') >/tmp/$USER/data
 *
 ./cns $(yad --columns=2 --title="Lorenz CNS (TSM)" --form --separator=" " --align=right \
    --field="Mode":CB --field="Separation" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="sigma" --field="rho" --field="beta (numerator)" --field="beta (denominator)" \
    -- 'step2!nosim' '1.0' './tsm-lorenz-static!./tsm-lorenz-dbg!./tsm-lorenz!./tsm-lorenz.py' \
    '6!3..64!3' '8!4..128!1' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' \
    '-15.8' '-18.48' '35.64' \
    '10' '28' '8' '3')
 *
 ./cns-scan $(yad --columns=2 --title="Lorenz CNS Scan (TSM)" --form --separator=" " --align=right \
    --field="Maxium Order":NUM --field="Separation" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":RO --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="sigma" --field="rho" --field="beta (numerator)" --field="beta (denominator)" \
    -- '32!2..64!1' '1.0' './tsm-lorenz-static!./tsm-lorenz-dbg!./tsm-lorenz!./tsm-lorenz.py' \
    '6!3..64!3' '_' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' \
    '-15.8' '-18.48' '35.64' \
    '10' '28' '8' '3') | tee /tmp/$USER/data
 *
 ./bifurcation-scan $(yad --columns=2 --title="Lorenz Bifurcation (TSM)" --form --separator=" " --align=right \
    --field="Lower Value" --field="Upper Value" --field="Skip Transient" --field="Model:CB" \
    --field="Display Places":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="sigma" --field="rho:RO" --field="beta (numerator)" --field="beta (denominator)" \
    -- '0.0' '50.0' '10' './tsm-lorenz-static!./tsm-lorenz-dbg!./tsm-lorenz' \
    '6!3..64!3' '4!4..128!1' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' \
    '-15.8' '-18.48' '35.64' \
    '10' '$p' '8' '3')
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "taylor-ode.h"

typedef struct { real sigma, rho, beta; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 12);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    real _;
    t_params(argv, argc, &p->sigma, &p->rho, &p->beta, &_);
    p->beta /= _;
    return p;
}

components ode (series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->sigma * (y[k] - x[k]),
        .y = p->rho * x[k] - y[k] - t_mul(x, z, k),
        .z = t_mul(x, y, k) - p->beta * z[k]
    };
}
