/*
 * Interface for performing high order Automatic Differentiation
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

#include "real.h"

/*
 * Initialize file scoped temporary storage
 */
void ad_init (int n);

/*
 * Set up a constant jet of value a, result stored in jet C
 */
series ad_const (series C, real a);

/*
 * Scales Taylor Series U by a factor a, result stored in jet S
 */
series ad_scale (series S, series U, real a);

/*
 * Calculates Taylor Series for the sum of U and V, result stored in jet P
 */
series ad_add (series P, series V, series U);

/*
 * Calculates Taylor Series for the difference of U and V, result stored in jet M
 */
series ad_sub (series M, series V, series U);

/*
 * Calculates Taylor Series for the absolute value of U, result stored in jet A
 */
series ad_abs (series A, series U);

/*
 * Calculates Taylor Series for U * V, result stored in jet P
 */
series ad_mul (series P, series V, series U);

/*
 * Calculates Taylor Series for U / V, result stored in jet Q
 */
series ad_div (series Q, series U, series V);

/*
 * Calculates Taylor Series for 1 / V, result stored in jet I
 */
series ad_inv (series I, series V);

/*
 * Calculates Taylor Series for the square of U, result stored in jet S
 */
series ad_sqr (series S, series U);

/*
 * Calculates Taylor Series for the square root of U, result stored in jet R
 */
series ad_sqrt (series R, series U);

/*
 * Calculates Taylor Series for the exponential of U, results stored in jet E
 */
series ad_exp (series E, series U);

/*
 * Calculates Taylor Series for the trigonometric sine and cosine of U, results stored in jets S and C
 */
void ad_sin_cos (series S, series C, series U, geometry g);

/*
 * Calculates Taylor Series for the trigonometric tangent and squared secant of U, results stored in jets T and S2
 */
void ad_tan_sec2 (series T, series S2, series U, geometry g);

/*
 * Calculates Taylor Series for U^a, results stored in jet P
 */
series ad_pwr (series P, series U, real a);

/*
 * Calculates Taylor Series for the natural logarithm of U, results stored in jet L
 */
series ad_ln (series L, series U);

void ad_asin (series AS, series DU_DF, series U, geometry G);

void ad_acos (series AC, series DU_DF, series U, geometry G);

void ad_atan (series AT, series DU_DF, series U, geometry G);
