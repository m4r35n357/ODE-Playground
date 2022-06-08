/*
 * Lorenz System
 *
 * Example: ./rk4-lorenz-dbg 15 1 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
 time -p ./cns step2 1 \
./rk4-lorenz-static $(yad --title="Lorenz CNS (RK4)" --form --separator=" " --align=right \
    --field="Display Precision":NUM --field="Plot Interval":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="sigma" --field="rho" --field="beta (numerator)" --field="beta (denominator)" \
    -- '6!3..64!3' '1!1..1000!1' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' "-15.8" "-18.48" "35.64" "10" "28" "8" "3")
 *
 time -p ./bifurcation-scan 0 50 10 \
./rk4-lorenz-static $(yad --title="Lorenz Bifurcation (RK4)" --form --separator=" " --align=right \
    --field="Display Precision":NUM --field="Plot Interval":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="sigma" --field="rho:RO" --field="beta (numerator)" --field="beta (denominator)" \
    -- '6!3..64!3' '1!1..1000!1' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' "-15.8" "-18.48" "35.64" "10" '$p' "8" "3")
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
