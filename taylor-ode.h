/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#pragma once
#include <time.h>
#include <stdbool.h>
#include "real.h"

/*
 * Client model data
 */
typedef struct Parameters model;

/*
 * Main floating point type
 */
typedef mpfr_t real;

/*
 * Type for Taylor Series coordinate jets
 */
typedef real *series;

/*
 * Combined x, y, z series
 */
typedef struct ts3 {
    series x, y, z;
} xyz;

/*
 * For returning x, y, z velocities from the model
 */
typedef struct triple_m {
    real x, y, z;
} triplet;

/*
 * Prints a line of data to stdout
 */
void _out_ (real x, real y, real z, real step_size, int step, clock_t since);

/*
 * Retrieves ODE parameters from the tail of the command (arguments 9 onwards)
 */
void tsm_get_p (char **argv, int count, ...);

/*
 * Creates a Taylor Series "jet" with the specified number of elements
 */
series tsm_jet (int size);

/*
 * Initialize constants
 */
void tsm_init (int display_precision);

/*
 * Safely and efficiently evaluates a polynomial of degree n, with the coefficients in S, and the variable in h
 */
real *horner (const series U, int order, real step_size);

/*
 *  Run TSM, send data to stdout
 */
void tsm (int order, real step_size, int steps, triplet *V, xyz *JETS, const model *p, clock_t since);

/*
 * Populate parameter data from command arguments
 */
model *tsm_init_p (int argc, char **argv, int order);

/*
 * Calculate kth components of the velocity jet V, using the ODE model together with the functions below as necessary.
 */
void ode (triplet *V, series X, series Y, series Z, const model *p, int k);

/*
 * Returns a pointer to kth element of the absolute value of U, no jet storage needed
 */
real *t_abs (const series U, int k);

/*
 * Returns a pointer to kth element of the product of U and V
 */
real *t_mul (const series U, const series V, const int k);

/*
 * Returns a pointer to kth element of the square of U
 */
real *t_sqr (const series U, int k);

/*
 * Calculates kth element of U / V, results stored in QUOT
 */
void t_div (series QUOT, const series U, const series V, int k);

/*
 * Calculates kth element of 1 / V, results stored in REC,
 */
void t_rec (series REC, const series V, int k);

/*
 * Calculates kth element of the square root of U, results stored in ROOT
 */
void t_sqrt (series ROOT, const series U, int k);

/*
 * Calculates kth element of the exponential of U, results stored in EXP
 */
void t_exp (series EXP, const series U, int k);

/*
 * Calculates kth elements of both sine and cosine of U, results stored in SIN and COS
 */
void t_sin_cos (series SIN, series COS, const series U, int k, bool trig);

/*
 * Calculates kth elements of both tan(h) and sec(h)^2 of U, results stored in TAN(H) and SEC(H)2
 */
void t_tan_sec2 (series TAN, series SEC2, const series U, int k, bool trig);

/*
 * Calculates kth element of the logarithm (inverse of EXP), results stored in U
 */
void t_ln (series U, const series EXP, int k);

/*
 * Calculates kth elements of arcsin/arsinh (inverse of SIN/SINH), results stored in U and COSH
 */
void t_asin_cos (series U, series COSH, const series SINH, int k, bool trig);

/*
 * Calculates kth elements of arccos/arcosh (inverse of COS/COSH), results stored in U and SINH
 */
void t_acos_sin (series U, series SINH, const series COSH, int k, bool trig);

/*
 * Calculates kth elements of arctan/artanh (inverse of TAN/TANH), results stored in U and SEC2
 */
void t_atan_sec2 (series U, series SEC2, const series TAN, int k, bool trig);

/*
 * Calculates kth element of P = U^a (where a is scalar), results stored in PWR
 */
void t_pwr (series PWR, series const U, real a, int k);
