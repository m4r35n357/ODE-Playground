/*
 * Thomas' cyclically symmetric attractor
 *
 * Example: ./tsm-thomas-dbg  9 237 16  0.1 30000  1 0 0  .185
 * 
 $(yad --columns=2 --title="Thomas Attractor (TSM)" --form --separator=" " --align=right \
    --field="Model:CB" \
    --field="Display Places":NUM --field="Precision (bits)":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="b" \
    -- './tsm-thomas-static!./tsm-thomas-dbg' \
    '6!3..64!3' '237!11..999!2' '8!2..256!1' '.1!0.001..0.1!0.001!3' '30000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '0.185') >/tmp/$USER/data
 *
 ./cns $(yad --columns=2 --title="Thomas CNS (TSM)" --form --separator=" " --align=right \
    --field="Mode":CB --field="Deviation" --field="Model:CB" \
    --field="Display Places":NUM --field="Precision (bits)":NUM --field="Order":NUM --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="b" \
    -- 'step2!nosim' '1.0' './tsm-thomas-static!./tsm-thomas-dbg' \
    '6!3..64!3' '237!11..999!2' '8!2..256!1' '.1!0.001..0.1!0.001!3' '30000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '0.185')
 *
 ./cns-scan $(yad --columns=2 --title="Thomas CNS Scan (TSM)" --form --separator=" " --align=right \
    --field="Minium Order":NUM --field="Maxium Order":NUM --field="Deviation" --field="Model:CB" \
    --field="Display Places":NUM --field="Precision (bits)":NUM --field="Order":RO --field="Step Size":NUM --field="Steps":NUM \
    --field="x0" --field="y0" --field="z0" \
    --field="b" \
    -- '2!2..256!1' '32!2..256!1' '1.0' './tsm-thomas-static!./tsm-thomas-dbg' \
    '6!3..64!3' '237!11..999!2' '_' '.1!0.001..0.1!0.001!3' '30000!1..1000000!1000' \
    '1.0' '0.0' '0.0' \
    '0.185') | tee /tmp/$USER/data
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

typedef struct { mpfr_t b; series sx, sy, sz, cx, cy, cz; } parameters;

void *get_p (int argc, char **argv, int n) {
    assert(argc == 10);
    parameters *p = malloc(sizeof (parameters));
    t_params(argv, argc, &p->b);
    p->sx = t_jet(n); p->cx = t_jet(n);
    p->sy = t_jet(n); p->cy = t_jet(n);
    p->sz = t_jet(n); p->cz = t_jet(n);
    return p;
}

void ode (components *vk, series x, series y, series z, void *params, int k) {
    parameters *p = (parameters *)params;
    //  x' = sin(y) - Bx
    mpfr_fms(vk->x, p->b, x[k], *t_sin_cos(p->sy, p->cy, y, k, TRIG).a, RND);
    mpfr_neg(vk->x, vk->x, RND);
    //  y' = sin(z) - By
    mpfr_fms(vk->y, p->b, y[k], *t_sin_cos(p->sz, p->cz, z, k, TRIG).a, RND);
    mpfr_neg(vk->y, vk->y, RND);
    //  z' = sin(x) - Bz
    mpfr_fms(vk->z, p->b, z[k], *t_sin_cos(p->sx, p->cx, x, k, TRIG).a, RND);
    mpfr_neg(vk->z, vk->z, RND);
}
