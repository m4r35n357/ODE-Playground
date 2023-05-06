/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

#include <time.h>
#include <stdbool.h>
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
controls *get_c_tsm (int argc, char **argv);

/*
 * Inital values for Taylor Series'
 */
series3 *initial_values (char **argv, int order);

/*
 * Retrieves ODE parameters from the tail of the command (arguments 8 onwards)
 */
void t_params (char **argv, int count, ...);

/*
 * Prints a line of data to stdout, with turning point markers
 */
void t_out (int dp, real x, real y, real z, real t, char *x_label, char *y_label, char *z_label, clock_t since);

/*
 * Creates a Taylor Series "jet" with the specified number of elements
 */
series t_jet (int size);

/*
 * Set up a constant jet of value a
 */
series t_const (int size, real a);

/*
 * Safely and efficiently evaluates a polynomial of degree n, with the coefficients in S, and the variable in h
 */
real t_horner (series S, int n, real h);

/*
 *  Dump data to stdout
 */
void tsm_stdout (int dp, controls *cont, series3 *jets, void *P, clock_t since);

/*
 * Generator implementation of TSM
 */
bool tsm_gen (controls *cont, series3 *jets, void *P);

/*
 * Obligatory client method signatures
 */

/*
 * Get a blob of parameter data from the model to be passed into ode()
 */
void *get_p (int argc, char **argv, int order);

/*
 * Calculate the kth components of the velocity jet V, using the coordinate jets and the parameter data,
 * together with the functions below as necessary.
 */
triplet ode (series X, series Y, series Z, void *P, int k);

/*
 * Basic Taylor Series functions
 */

/*
 * Returns kth element of the absolute value of U, no user-supplied jet storage needed
 */
real t_abs (series U, int k);

/*
 * Taylor Series recurrence relationships
 */

/*
 * Returns kth element of the product of U and V, no user-supplied jet storage needed
 */
real t_mul (series U, series V, int k);

/*
 * Returns kth element of the square of U, no user-supplied jet storage needed
 */
real t_sqr (series U, int k);

/*
 * Returns kth element of U / V, results stored in user-supplied jet Q, DOMAIN RESTRICTION v[0] != 0.0
 */
real t_div (series Q, series U, series V, int k);

/*
 * Returns kth element of the square root of U, results stored in user-supplied jet R, DOMAIN RESTRICTION U[0] > 0.0
 */
real t_sqrt (series R, series U, int k);

/*
 * Returns kth element of the exponential of U, results stored in user-supplied jet E
 */
real t_exp (series E, series U, int k);

/*
 * Returns kth elements of both sine and cosine of U, results stored in user-supplied jets S and C
 */
pair t_sin_cos (series S, series C, series U, int k, bool trig);

/*
 * Returns kth elements of both tangent and squared secant of U, results stored in user-supplied jets T and S2
 */
pair t_tan_sec2 (series T, series S2, series U, int k, bool trig);

/*
 * Returns kth element of P = U^a (where a is scalar), results stored in user-supplied jet P, DOMAIN RESTRICTION U[0] > 0.0
 */
real t_pwr (series P, series U, real a, int k);

/*
 * Returns kth element of the natural logarithm of U, results stored in user-supplied jet L, DOMAIN RESTRICTION U[0] > 0.0
 */
real t_ln (series L, series U, int k);

/*
 * Returns kth elements of arcsin(h) of U and 1 / DF_DU, results stored in user-supplied jets As and DU_DF
 */
pair t_asin (series AS, series DU_DF, series U, int k, bool trig);

/*
 * Returns kth elements of arccos(h) of U and 1 / DF_DU, results stored in user-supplied jets As and DU_DF
 */
pair t_acos (series AC, series G, series U, int k, bool trig);

/*
 * Returns kth elements of arctan(h) of U and 1 / DF_DU, results stored in user-supplied jets As and DU_DF
 */
pair t_atan (series AT, series G, series U, int k, bool trig);
