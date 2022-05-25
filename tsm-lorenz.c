/*
 * Lorenz System
 *
 * Example: ./tsm-lorenz-dbg 9 237 16 .01 10000 -15.8 -17.48 35.64 10 28 8 3
 *
./cns step2 1 ./tsm-lorenz-static $(yad --title="Lorenz Attractor (TSM)" --form --separator=" " \
    --field="Display Precision":NUM \
    --field="Precision in Bits":NUM \
    --field="Order":NUM \
    --field="Step Size":NUM \
    --field="Steps":NUM \
    --field="x0" \
    --field="y0" \
    --field="z0" \
    --field="Sigma" \
    --field="Rho" \
    --field="Beta numerator" \
    --field="Beta denominator" \
    -- '6!3..64!3' '237!11..999!2' '8!4..256!1' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' "-15.8" "-18.48" "35.64" "10" "28" "8" "3")
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t sigma, rho, beta, _; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 13);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->sigma, &p->rho, &p->beta, &p->_);
    mpfr_div(p->beta, p->beta, p->_, RND);
    return p;
}

void ode (components *vk, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = S(y - x)
    mpfr_fmms(vk->x, p->sigma, y[k], p->sigma, x[k], RND);
    //  y' = x(R - z) - y
    mpfr_fms(vk->y, x[k], p->rho, *t_mul(x, z, k), RND);
    mpfr_sub(vk->y, vk->y, y[k], RND);
    //  z' = xy - Bz
    mpfr_fms(vk->z, p->beta, z[k], *t_mul(x, y, k), RND);
    mpfr_neg(vk->z, vk->z, RND);
}
