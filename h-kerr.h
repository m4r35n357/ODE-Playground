/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

#include "dual.h"

typedef struct Kerr {
    real mu2;  // central mass & particle mass (squared)
    real E, L, Q, K;  // constants of motion
    real a, a2, L2, aL, aE, a2xmu2_E2;  // global constants
    real step_size, tau, q_t, q_r, q_theta, q_phi, v_t, v_r, v_theta, v_phi;  // proper time, coordinates & velocities
    dual ra2, D, sth2, R, TH;  // global variables & potentials
    components *coordinates;
    real rmin, rmax, thmax, epsilon;  // generator constraints and precision
    float horizon;
} kerr;

/*
 * Get a blob of model data from the command to be passed into solve()
 */
kerr *get_p_kerr (int argc, char **argv, real step);

real elevation_to_colatitude (real elevation);

real sigma (kerr *bh);

pair gamma_v (kerr *bh, real sigma);
