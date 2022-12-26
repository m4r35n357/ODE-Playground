/*
 * N-body problem using Hamilton's equations
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "symplectic.h"
#include "h-nbody.h"

static void plot (int dp, void *n_body, real t) {
    nbody *nb = (nbody *)n_body;
    body *b = nb->bodies;
    printf("%.6Le %+*Le", t, dp, error(h(nb) - nb->h0));
    for (int i = 0; i < nb->n; i++) {
        printf("  %+*Le %+*Le %+*Le %+*Le %+*Le %+*Le",
               dp, b[i].x, dp, b[i].y, dp, b[i].z, dp, b[i].px, dp, b[i].py, dp, b[i].pz);
    }
    printf("\n");
}

int main (int argc, char** argv) {
    solve(argv, get_c_symp(argc, argv), get_p_nbody(argc, argv, (argc - 6) / 7), plot);
    return 0 ;
}
