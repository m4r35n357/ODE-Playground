
#include "real.h"

typedef struct Body {
    real m;  // mass
    real q_x, q_y, q_z, p_x, p_y, p_z;  // coordinates & momenta
} body;

typedef struct Nbody {
    int n;
    body *bodies;
    real g, h0;
} nbody;

typedef struct Components {
    real x;
    real y;
    real z;
} components;

components cog (nbody *nb);
