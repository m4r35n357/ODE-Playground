/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

#include "dual.h"

struct Parameters {
    real mu2;  // central mass & particle mass (squared)
    real E, L, Q, K;  // constants of motion
    real a, a2, L2, aL, aE, a2xmu2_E2;  // global constants
    real step_size, tau, q_t, q_r, q_th, q_ph, v_t, v_r, v_th, v_ph;  // proper time, coordinates & velocities
    dual ra2, D, sth2, R, TH;  // global variables & potentials
    triplet *coordinates;
    real rmin, rmax, thmax, epsilon;  // generator constraints and precision
    float horizon;
};

/*
 * Get model data from the command
 */
model *kerr_get_p (int argc, char **argv, real step);

real elevation_to_colatitude (real elevation);

real sigma (model *bh);

pair gamma_v (model *bh, real sigma);
