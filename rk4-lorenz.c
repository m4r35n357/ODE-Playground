/*
 * Lorenz System
 *
 * Example: ./rk4-lorenz-dbg  6 1 .01 10000  -15.8 -17.48 35.64  10 28 8 3
 *
 $(yad --columns=2 --title="Lorenz Attractor (RK4)" --form --separator=" " --align=right \
    --field="Model:CB" \
    --field="Display Precision":NUM --field="Plot Interval":NUM --field="Step Size" --field="Steps" \
    --field="x0" --field="y0" --field="z0" \
    --field="sigma" --field="rho" --field="beta (numerator)" --field="beta (denominator)" \
    -- './rk4-lorenz-static!./rk4-lorenz-dbg!./rk4-lorenz' \
    '6!3..64!3' '1!1..1000!1' '.01' '10000' \
    '-15.8' '-18.48' '35.64' \
    '10' '28' '8' '3') >/tmp/$USER/data
 *
 ./cns $(yad --columns=2 --title="Lorenz CNS (RK4)" --form --separator=" " --align=right \
    --field="Mode":CB --field="Separation" --field="Model:CB" \
    --field="Display Precision":NUM --field="Plot Interval":NUM --field="Step Size" --field="Steps" \
    --field="x0" --field="y0" --field="z0" \
    --field="sigma" --field="rho" --field="beta (numerator)" --field="beta (denominator)" \
    -- 'step2!nosim' '1' './rk4-lorenz-static!./rk4-lorenz-dbg!./rk4-lorenz' \
    '6!3..64!3' '1!1..1000!1' '.01' '10000' \
    '-15.8' '-18.48' '35.64' \
    '10' '28' '8' '3')
 *
 ./bifurcation-scan $(yad --columns=2 --title="Lorenz Bifurcation (RK4)" --form --separator=" " --align=right \
    --field="Lower Value" --field="Upper Value" --field="Skip Transient" --field="Model:CB" \
    --field="Display Precision":NUM --field="Plot Interval":NUM --field="Step Size" --field="Steps" \
    --field="x0" --field="y0" --field="z0" \
    --field="sigma" --field="rho:RO" --field="beta (numerator)" --field="beta (denominator)" \
    -- '0.0' '50.0' '10' './rk4-lorenz-static!./rk4-lorenz-dbg!./rk4-lorenz' \
    '6!3..64!3' '1!1..1000!1' '.01' '10000' \
    '-15.8' '-18.48' '35.64' \
    '10' '$p' '8' '3')
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include "runge-kutta-ode.h"

typedef struct { real sigma, rho, beta; } parameters;

void *get_p (int argc, char **argv) {
    assert(argc == 12);
    parameters *p = malloc(sizeof (parameters));
    real _;
    t_params(argv, argc, &p->sigma, &p->rho, &p->beta, &_);
    p->beta /= _;
    return p;
}

components ode (real x, real y, real z, void *params) {
    parameters *p = (parameters *)params;
    return (components) {
        .x = p->sigma * (y - x),
        .y = p->rho * x - y - x * z,
        .z = x * y - p->beta * z
    };
}
