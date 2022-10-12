
#include "dual.h"

typedef struct Parameters {
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
} parameters;

point to_xyz (parameters *bh);

real r (parameters *bh);
real theta (parameters *bh);
real phi (parameters *bh);

real sigma (parameters *bh);

pair gamma_v (parameters *bh, real sigma);
