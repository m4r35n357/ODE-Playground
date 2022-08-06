
#include "dual.h"

typedef struct Parameters {
    real mu2;  // central mass & particle mass (squared)
    real E, L, Q, K;  // constants of motion
    real a, a2, L2, aL, aE, a2xmu2_E2;  // global constants
    real q_t, q_r, q_theta, q_phi, p_t, p_r, p_theta, p_phi;  // coordinates & velocities
    dual ra2, delta, sth2, R, THETA;  // global variables & potentials
    components *coordinates;
    rgb colour, *track;
    float ball_size, view_radius, view_latitude, view_longitude, horizon;
    int max_points, oldest, newest;
    _Bool buffers_full;
} parameters;

rgb to_xyz (parameters *bh);
real r (parameters *bh);
real theta (parameters *bh);
real phi (parameters *bh);
pair gamma (parameters *bh);

void plot_path (int dp, void *params, real t);
void plot_view (int dp, void *params, real t);
void plot_raw (int dp, void *params, real time);
