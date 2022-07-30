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
