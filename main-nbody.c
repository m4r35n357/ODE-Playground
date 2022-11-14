/*
 * N-Body simulator
 *
 * Example:  ./h-nbody-std 6 8 0.010 10000 _ .05 100.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 0.0 4.5 0.4 -0.2 0.0 1.8 3.0 -6.0 0.0 -0.4 0.0 -2.0 1.0 5.0 3.0 0.0 -0.2 0.0 5.8 -0.2 4.0 0.0 -4.0 0.1 -3.6 0.0 0.2 3.0 -4.0 0.0 -0.1 0.0 -0.2 -2.6 3.0 8.0 0.0 -0.3 0.0 2.0 -0.2 4.0 0.0 4.0 -0.2 -4.8 0.0 -0.2
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include "symplectic.h"
#include "h-nbody.h"

static void plot (int dp, void *n_body, real t) {
    nbody *nb = (nbody *)n_body;
    body *b = nb->bodies;
    printf("%.6Le", t);
    for (int i = 0; i < nb->n; i++) {
        printf("  %+*Le %+*Le %+*Le %+*Le %+*Le %+*Le",
               dp, b[i].x, dp, b[i].y, dp, b[i].z, dp, b[i].px, dp, b[i].py, dp, b[i].pz);
    }
    printf("\n");
}

int main (int argc, char** argv) {
    solve(argv, get_c_symp(argv), get_p_nbody(argc, argv, (argc - 7) / 7), plot);
    return(0);
}
