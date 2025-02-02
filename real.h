/*
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#pragma once
#include <stdbool.h>

/*
 * Number base for integer conversions
 */
#define BASE 10

/*
 * Colours
 */
#define GRY "\x1B[1;30m"
#define RED "\x1B[1;31m"
#define GRN "\x1B[1;32m"
#define YLW "\x1B[1;33m"
#define BLU "\x1B[1;34m"
#define MGT "\x1B[0;35m"
#define CYN "\x1B[0;36m"
#define WHT "\x1B[1;37m"
#define NRM "\x1B[0m"

/*
 * "Low-noise" squaring for arguments with no side-effects
 */
#define SQR(x) ((x) * (x))

/*
 * Report program arguments in colour
 */
#define PRINT_ARGS(argc, argv) do { \
    fprintf(stderr, "%sargc %s%2d%s, argv [ %s", GRY, NRM, (argc), GRY, CYN); \
    for (int i = 0; i < (argc); i++) fprintf(stderr, "%s ", (argv)[i]); \
    fprintf(stderr, "%s]%s\n", GRY, NRM); \
} while (0)

/*
 * Unavoidable "assert", in colour
 */
#define CHECK(x) do { \
    if(!(x)) { \
        fprintf(stderr, "%sFAIL %s%s %s%s%s %s%s:%s%i\n", RED, WHT, #x, GRY, __func__, NRM, __FILE__, GRY, NRM, __LINE__); \
        abort(); \
    } \
} while (0)

/*
 * Main floating point type
 */
typedef long double real;

/*
 * Client model data
 */
typedef struct Parameters model;

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
    bool looping;         // generators only
    int order, step, steps, dp;
    real h;
} controls;
