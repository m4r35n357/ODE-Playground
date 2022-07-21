
#include "real.h"

typedef struct Components {
    real x, y, z;
} components;

typedef struct Line {
    int newest;
    components *buffer;
} line;

typedef struct Rgb {
    float r, g, b;
} rgb;

typedef struct Body {
    real m;  // mass
    real q_x, q_y, q_z, p_x, p_y, p_z;  // coordinates & momenta
    rgb colour;
    line *track;
} body;

typedef struct Nbody {
    int n;
    body *bodies;
    real g, h0;
    components centre;
    float ball_scale, view_radius, view_latitude, view_longitude;
} nbody;

void cog (nbody *nb);

real h (nbody *nb);
