
#include "dual.h"

typedef struct Kerr {
    real mu2;  // central mass & particle mass (squared)
    real E, L, Q, K;  // constants of motion
    real a, a2, L2, aL, aE, a2xmu2_E2;  // global constants
    real step, tau, q_t, q_r, q_theta, q_phi, p_t, p_r, p_theta, p_phi;  // proper time, coordinates & velocities
    dual ra2, delta, sth2, R, THETA;  // global variables & potentials
    components *coordinates;
    struct triple_f colour, *track;
    float ball_scale, view_radius, view_latitude, view_longitude, horizon;
    int max_points, oldest, newest;
    _Bool buffers_full;
    real rmin, rmax, thmax, epsilon;  // generator constraints and precision
} kerr;

/*
 * Get a blob of model data from the command to be passed into solve()
 */
kerr *get_p_kerr (int argc, char **argv);

real elevation_to_colatitude (real elevation);

real sigma (kerr *bh);

pair gamma_v (kerr *bh, real sigma);
