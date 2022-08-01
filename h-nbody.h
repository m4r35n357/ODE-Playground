
#include "real.h"

typedef struct Body {
    real m, x, y, z, px, py, pz;
    float r;
    rgb colour;
    rgb *track;
} body;

typedef struct Nbody {
    int n;
    body *bodies;
    real g, h0, h;
    components centre;
    float ball_scale, view_radius, view_latitude, view_longitude;
    int max_points, oldest, newest;
    _Bool buffers_full;
} nbody;

/*
 * Calculate the centre of mass of the system
 */
void cog (nbody *nb);

/*
 * Hamiltonian
 */
real h (nbody *nb);
