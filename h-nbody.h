
#include "real.h"

typedef struct Body {
    real m, x, y, z, px, py, pz;
    float r;
    struct triple_f colour, *track;
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
 * Get a blob of model data from the command to be passed into solve()
 */
nbody *get_p_nbody (int argc, char **argv, int va_begin);

/*
 * Calculate the centre of mass of the system
 */
void cog (nbody *nb);

/*
 * Hamiltonian
 */
real h (nbody *nb);
