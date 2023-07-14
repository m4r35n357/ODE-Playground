/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

typedef long double real;

/*
 * Client model data
 */
typedef struct Parameters parameters;

/*
 * For returning combined values from trig/hyp functions
 */
typedef struct pair_l {
    real a, b;
} pair;

/*
 * Triple of long doubles, by coordinate
 */
typedef struct triple_l {
    real x, y, z;
} triplet;

/*
 * Integrator parameters
 */
typedef struct Controls {
    int order, step, steps;
    real step_size;
    pair r2, r4, r6, r8;  // symplectic only
} controls;

/*
 * Number base for integer conversions
 */
#define BASE 10

/*
 * Report program arguments in colour
 */
#define PRINT_ARGS(argc, argv) \
    fprintf(stderr, "argc: \x1B[1;37m%d\x1B[0;37m, argv: [ \x1B[0;36m", argc); \
    for (int i = 0; i < argc; i++) { \
        fprintf(stderr, "%s ", argv[i]); \
    } \
    fprintf(stderr, "\x1B[0;37m]\n");

/*
 * Unavoidable assert(), in colour
 */
#define CHECK(x) \
    do { \
        if(!(x)) { \
            fprintf(stderr, \
                "\x1B[1;31mFAIL\x1B[0;37m \x1B[1;37m%s\x1B[0;37m \x1B[0;35m%s\x1B[0;37m:\x1B[0;35m%i\x1B[0;37m\n", \
                #x, __FILE__, __LINE__); \
            exit(1); \
        } \
    } while (0)
