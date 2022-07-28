/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

typedef long double real;

static const int BASE = 10;

static const real MY_PI = 3.1415926535897932384626433832795029L;

/*
 * General x, y, z components
 */
typedef struct Components {
    real x, y, z;
} components;

/*
 * Colour r, g, b values
 */
typedef struct Rgb {
    float r, g, b;
} rgb;

/*
 * Particle/Body tracks
 */
typedef struct Line {
	int newest;
    components *buffer;
} line;

typedef struct Particle {
    components *coordinates;
    rgb colour;
    line *track;
    float size, view_radius, view_latitude, view_longitude;
} particle;

typedef struct Body {
    real m;
    real x, y, z, px, py, pz;
    rgb colour;
    line *track;
} body;

typedef struct Nbody {
    int n;
    body *bodies;
    real g, h0, h;
    components centre;
    float ball_scale, view_radius, view_latitude, view_longitude;
    int max_points, oldest, newest;
    _Bool full;
} nbody;

/*
 * Integrator control parameters
 */
typedef struct Weights {
    real fwd, rev;
} weights;

typedef struct Controls {
    int order, step, steps;
    real step_size;
    weights r1, r2, r3, r4;
} controls;

/*
 * Retrieves integrator control parameters
 */
controls *get_c (char **argv);
