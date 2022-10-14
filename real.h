/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

typedef long double real;

static const int BASE = 10;

static const real MY_PI = 3.1415926535897932384626433832795029L;

static const real RAD_TO_DEG = 180.0L / MY_PI;

/*
 * For returning combined values
 */
typedef struct pair_l {
    real a, b;
} pair;

/*
 * General x, y, z components
 */
typedef struct triple_l {
    real x, y, z;
} components;

/*
 * Colour r, g, b values
 */
typedef struct triple_f {
    float a, b, c;
} rgb, point;

/*
 * Integrator control parameters
 */
typedef struct Weights {
    real fwd, rev;
} weights;

typedef struct Controls {
    int order, step, steps, generating;
    real step_size;
    weights r1, r2, r3, r4;
} controls;
