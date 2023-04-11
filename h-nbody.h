/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

#include "real.h"

typedef struct Body {
    real m, x, y, z, px, py, pz;
    float r;
} body;

typedef struct Nbody {
    int n;
    body *bodies;
    real g, h0;
} nbody;

/*
 * Get a blob of model data from the command to be passed into solve()
 */
nbody *get_p_nbody (int argc, char **argv);

/*
 * Calculate the centre of mass of the system
 */
void reset_cog (nbody *nb);

/*
 * Hamiltonian
 */
real hamiltonian (nbody *nb);
