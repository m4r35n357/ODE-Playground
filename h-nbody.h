/*
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#pragma once
#include "real.h"

typedef struct Body {
    real m, x, y, z, px, py, pz;
    float r;
} body;

struct Parameters {
    int n;
    body *bodies;
    real G, h0;
};

/*
 * Get a blob of model data from the command to be passed into solve()
 */
model *get_p_nbody (int argc, char **argv);

/*
 * Calculate the centre of mass of the system
 */
void reset_cog (model *nb);

/*
 * Hamiltonian
 */
real H (model *nb);
