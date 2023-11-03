/*
 * N-body problem using Hamilton's equations
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "symplectic.h"
#include "h-nbody.h"

static void plot (int dp, model *nb, real t) {
    body *b = nb->bodies;
    reset_cog(nb);
    printf("%.6Le %+.*Le", t, dp, error(H(nb) - nb->h0));
    for (int i = 0; i < nb->n; i++) {
        printf("  %+.*Le %+.*Le %+.*Le %+.*Le %+.*Le %+.*Le",
                dp, b[i].x, dp, b[i].y, dp, b[i].z, dp, b[i].px, dp, b[i].py, dp, b[i].pz);
    }
    printf("\n");
}

int main (int argc, char **argv) {
    solve(symp_get_c(argc, argv), get_p_nbody(argc, argv), plot);
    return 0 ;
}
