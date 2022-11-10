/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

extern const int BASE;

typedef long double real;

/*
 * For returning combined values
 */
typedef struct pair_l {
    real a, b;
} pair;

/*
 * Triple of long doubles
 */
typedef struct triple_l {
    real x, y, z;
} components;

/*
 * Triple of floats
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
    int order, step, steps;
    real step_size;
    weights r2, r4, r6, r8;
} controls;
