/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
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

#define PRINT_ARGS(argc, argv) \
    fprintf(stderr, "[ "); for (int i = 0; i < argc; i++) fprintf(stderr, "%s ", argv[i]); fprintf(stderr, "]\n");
