
#include "real.h"

typedef struct Components {
    real x;
    real y;
    real z;
} components;

typedef struct Body {
    real m;  // mass
    real q_x, q_y, q_z, p_x, p_y, p_z;  // coordinates & momenta
    components colour;
} body;

typedef struct Nbody {
    int n;
    body *bodies;
    real g, h0;
    components centre;
    float radius, view_latitude, view_longitude;
} nbody;

void cog (nbody *nb);
