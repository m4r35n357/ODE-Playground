/*
 * Interface for solving systems of Ordinary Differential Equations using the Taylor Series Method (TSM)
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
typedef struct ts3 {
    series x, y, z;
} xyz;

/*
 * Retrieves control parameters from the first four command arguments
 */
controls *tsm_get_c (int argc, char **argv);

/*
 * Retrieves initial X, Y, Z values from the next three command arguments and populates their Taylor Series'
 */
xyz *tsm_init (char **argv, int order);

/*
 * Retrieves ODE model parameters from the tail of the command (argument 8 onwards)
 */
void tsm_get_p (char **argv, int count, ...);

/*
 * Creates a Taylor Series with the specified number of elements
 */
series tsm_jet (int size);

/*
 * Safely and efficiently evaluates a polynomial of degree n, with the coefficients in S, and the variable in h
 */
real horner (series U, int order, real h);

/*
 *  Run TSM, send data to stdout
 */
void tsm_out (controls *c, xyz *jets, const model *p, clock_t since);

/*
 * Generator (step-wise) implementation of TSM
 */
bool tsm_gen (controls *c, xyz *jets, const model *p);

/*
 * Obligatory client method signatures
 */

/*
 * Populate parameter data from command arguments
 */
model *tsm_init_p (int argc, char **argv, int order);

/*
 * Calculate kth components of the velocity V, using the ODE model together with the functions below as necessary.
 */
triplet ode (series X, series Y, series Z, const model *p, int k);

/*
 * Basic Taylor Series functions
 */

/*
 * Returns kth element of the absolute value of U, no storage needed
 */
real t_abs (const series U, int k);

/*
 * Taylor Series recurrence relationships
 */

/*
 * Returns kth element of the product of U and W, no storage needed
 */
real t_mul (const series U, const series W, int k);

/*
 * Returns kth element of the square of U, no storage needed
 */
real t_sqr (const series U, int k);

/*
 * Returns kth element of U / W, results stored in QUOT
 */
real t_div (series QUOT, const series U, const series W, int k);

/*
 * Returns kth element of the square root of U, results stored in ROOT
 */
real t_sqrt (series ROOT, const series U, int k);

/*
 * Returns kth element of the exponential of U, results stored in EXP
 */
real t_exp (series EXP, const series U, int k);

/*
 * Returns kth elements of both SIN/SINH and COS/COSH of U, results stored in SIN/SINH and COS/COSH
 */
pair t_sin_cos (series SIN, series COS, const series U, int k, bool trig);

/*
 * Returns kth elements of both TAN/TANH and SEC2/SECH2 of U, results stored in TAN/TANH and SEC2/SECH2
 */
pair t_tan_sec2 (series TAN, series SEC2, const series U, int k, bool trig);

/*
 * Returns kth element of the logarithm (inverse of EXP), results stored in U
 */
real t_ln (series U, const series EXP, int k);

/*
 * Returns kth elements of arcsin/arsinh (inverse of SIN/SINH) with COS/COSH, results stored in U and COS/COSH
 */
pair t_asin_cos (series U, series COS, const series SIN, int k, bool trig);

/*
 * Returns kth elements of arccos/arcosh (inverse of COS/COSH) with SIN/SINH, results stored in U and SIN/SINH
 */
pair t_acos_sin (series U, series SIN, const series COS, int k, bool trig);

/*
 * Returns kth elements of arctan/artanh (inverse of TAN/TANH) with SEC2/SECH2, results stored in U and SEC2/SECH2
 */
pair t_atan_sec2 (series U, series SEC2, const series TAN, int k, bool trig);

/*
 * Returns kth element of P = U^a (where a is scalar), results stored in PWR
 */
real t_pwr (series PWR, series U, real a, int k);
