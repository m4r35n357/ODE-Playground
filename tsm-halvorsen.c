/*
 * Halvorsen Cyclic Attractor
 *
 * Example: ./tsm-halvorsen-dbg  6 237 16  .01 10000  1 0 0  1.4 4
 *
 $(yad --columns=2 --title="Halvorsen Attractor (TSM)" --form --separator=" " --align=right \
    --field="Model:CB" \
    --field="Display Places":NUM --field="Precision (bits)":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="a" --field="b" \
    -- './tsm-halvorsen-static!./tsm-halvorsen-dbg' \
    '6!3..64!3' '237!11..999!2' '16!2..256!1' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '1.4' '4') >/tmp/$USER/data
 *
 ./cns $(yad --columns=2 --title="Halvorsen CNS (TSM)" --form --separator=" " --align=right \
    --field="Mode":CB --field="Deviation" --field="Model:CB" \
    --field="Display Places":NUM --field="Precision":CB --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="a" --field="b" \
    -- 'step2!nosim' '1.0' './tsm-halvorsen-static!./tsm-halvorsen-dbg' \
    '6!3..64!3' 'octuple!quadruple!extended!double!single' '16!2..256!1' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '1.4' '4')
 *
 ./cns-scan $(yad --columns=2 --title="Halvorsen CNS Scan (TSM)" --form --separator=" " --align=right \
    --field="Minium Order":NUM --field="Maxium Order":NUM --field="Deviation" --field="Model:CB" \
    --field="Display Places":NUM --field="Precision":CB --field="Order":RO --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="a" --field="b" \
    -- '2!2..256!1' '32!2..256!1' '1.0' './tsm-halvorsen-static!./tsm-halvorsen-dbg' \
    '6!3..64!3' 'octuple!quadruple!extended!double!single' '_' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '1.4' '4') | tee /tmp/$USER/data
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t a, b; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 11);
    (void)n;
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->a, &p->b);
    return p;
}

void ode (components *vk, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = - Ax - By - Bz - y^2
    mpfr_fmma(vk->x, p->b, y[k], p->b, z[k], RND);
    mpfr_fma(vk->x, p->a, x[k], vk->x, RND);
    mpfr_add(vk->x, *t_sqr(y, k), vk->x, RND);
    mpfr_neg(vk->x, vk->x, RND);
    //  y' = - Ay - Bz - Bx - z^2
    mpfr_fmma(vk->y, p->b, z[k], p->b, x[k], RND);
    mpfr_fma(vk->y, p->a, y[k], vk->y, RND);
    mpfr_add(vk->y, *t_sqr(z, k), vk->y, RND);
    mpfr_neg(vk->y, vk->y, RND);
    //  z' = - Az - Bx - By - x^2
    mpfr_fmma(vk->z, p->b, x[k], p->b, y[k], RND);
    mpfr_fma(vk->z, p->a, z[k], vk->z, RND);
    mpfr_add(vk->z, *t_sqr(x, k), vk->z, RND);
    mpfr_neg(vk->z, vk->z, RND);
}
