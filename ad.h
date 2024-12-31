/*
 * Interface for performing high order Automatic Differentiation
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

/*
 * Initialize file scoped temporary storage
 */
void ad_init (int n);


/*
 * Scales Taylor Series U by a factor a, result stored in jet SCL
 */
series ad_scale (series SCL, series U, real a);

/*
 * Calculates Taylor Series for the sum of U and V, result stored in jet SUM
 */
series ad_add (series SUM, series V, series U);

/*
 * Calculates Taylor Series for the difference of U and V, result stored in jet DIFF
 */
series ad_sub (series DIFF, series V, series U);

/*
 * Calculates Taylor Series for the absolute value of U, result stored in jet ABS
 */
series ad_abs (series ABS, series U);

/*
 * Calculates Taylor Series for U * V, result stored in jet PROD
 */
series ad_mul (series PROD, series V, series U);

/*
 * Calculates Taylor Series for U / V, result stored in jet QUOT
 */
series ad_div (series QUOT, series U, series V);
series ad_rec (series r, series v);

/*
 * Calculates Taylor Series for the square of U, result stored in jet SQR
 */
series ad_sqr (series SQR, series U);

/*
 * Calculates Taylor Series for the square root of U, result stored in jet ROOT
 */
series ad_sqrt (series ROOT, series U);

/*
 * Calculates Taylor Series for the exponential of U, results stored in jet EXP
 */
series ad_exp (series EXP, series U);

/*
 * Calculates Taylor Series for the trigonometric sine and cosine of U, results stored in jets SIN and COS
 */
pair ad_sin_cos (series SIN, series COS, series U, bool trig);

/*
 * Calculates Taylor Series for the trigonometric tangent and squared secant of U, results stored in jets TAN and SEC2
 */
pair ad_tan_sec2 (series TAN, series SEC2, series U, bool trig);

/*
 * Calculates Taylor Series for the natural logarithm (inverse of EXP), results stored in jet U
 */
series ad_ln (series U, series EXP);

/*
 * Returns kth elements of arcsin/arsinh (inverse of SIN/SINH), results stored in jets U and COSH
 */
void ad_asin (series U, series COSH, series SINH, bool trig);

/*
 * Returns kth elements of arccos/arcosh (inverse of COS/COSH), results stored in jets U and SINH
 */
void ad_acos (series U, series SINH, series COSH, bool trig);

/*
 * Returns kth elements of arctan/artanh (inverse of TAN/TANH), results stored in jets U and SEC2
 */
void ad_atan (series U, series SEC2, series TAN, bool trig);

/*
 * Calculates Taylor Series for U^a, results stored in jet P
 */
series ad_pwr (series P, series U, real a);
