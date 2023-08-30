/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

#include <time.h>
#include "real.h"

/*
 * Type for Taylor Series coordinate jets
 */
typedef real *series;

/*
 * Combined x, y, z series
 */
typedef struct triple_s {
    series x, y, z;
} series3;

/*
 * Retrieves integrator control parameters
 */
controls *tsm_get_c (int argc, char **argv);

/*
 * Initial values for Taylor Series'
 */
series3 *tsm_init_xyz (char **argv, int order);

/*
 * Retrieves ODE parameters from the tail of the command (arguments 8 onwards)
 */
void tsm_get_p (char **argv, int count, ...);

/*
 * Creates a Taylor Series "jet" with the specified number of elements
 */
series tsm_var (int size);

/*
 * Set up a constant jet of value a
 */
series tsm_const (int size, real a);

/*
 * Safely and efficiently evaluates a polynomial of degree n, with the coefficients in S, and the variable in h
 */
real horner (series S, int n, real h);

/*
 *  Run TSM, send data to stdout
 */
void tsm_stdout (int dp, controls *cont, series3 *jets, parameters *P, clock_t since);

/*
 * Generator (step-wise) implementation of TSM
 */
bool tsm_gen (controls *cont, series3 *jets, parameters *P);

/*
 * Obligatory client method signatures
 */

/*
 * Populate parameter data from command arguments
 */
parameters *tsm_init_p (int argc, char **argv, int order);

/*
 * Calculate kth components of the velocity jet V, using the ODE model together with the functions below as necessary.
 */
triplet ode (series X, series Y, series Z, parameters *P, int k);

/*
 * Basic Taylor Series functions
 */

/*
 * Returns kth element of the absolute value of U, no jet storage needed
 */
real t_abs (series U, int k);

/*
 * Taylor Series recurrence relationships
 */

/*
 * Returns kth element of the product of U and W, no jet storage needed
 */
real t_mul (series U, series W, int k);

/*
 * Returns kth element of the square of U, no jet storage needed
 */
real t_sqr (series U, int k);

/*
 * Returns kth element of U / W, results stored in jet QUOT
 */
real t_div (series QUOT, series U, series W, int k);

/*
 * Returns kth element of the square root of U, results stored in jet ROOT
 */
real t_sqrt (series ROOT, series U, int k);

/*
 * Returns kth element of the exponential of U, results stored in jet EXP
 */
real t_exp (series EXP, series U, int k);

/*
 * Returns kth elements of both sine and cosine of U, results stored in jets SIN and COS
 */
pair t_sin_cos (series SIN, series COS, series U, int k, bool trig);

/*
 * Returns kth elements of both tangent and squared secant of U, results stored in jets TAN and SEC2
 */
pair t_tan_sec2 (series TAN, series SEC2, series U, int k, bool trig);

/*
 * Returns kth element of the logarithm (inverse of EXP), results stored in jet U
 */
real t_ln (series U, series EXP, int k);

/*
 * Returns kth elements of arcsin/arsinh (inverse of SIN/SINH), results stored in jets U and COSH
 */
pair t_asin (series U, series COSH, series SINH, int k, bool trig);

/*
 * Returns kth elements of arccos/arcosh (inverse of COS/COSH), results stored in jets U and SINH
 */
pair t_acos (series U, series SINH, series COSH, int k, bool trig);

/*
 * Returns kth elements of arctan/artanh (inverse of TAN/TANH), results stored in jets U and SEC2
 */
pair t_atan (series U, series SEC2, series TAN, int k, bool trig);

/*
 * Returns kth element of P = U^a (where a is scalar), results stored in jet PWR
 */
real t_pwr (series PWR, series U, real a, int k);
