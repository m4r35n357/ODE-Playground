
/*
 * Interface for performing high order Automatic Differentiation
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

/*
 * Whether a jet is intended for automatic diffentiation, or not
 */
typedef enum {VARIABLE, CONSTANT} ad_status;

/*
 * Solver method
 */
typedef enum {NONE=0, NEWTON=2} solver;

/*
 * Solver mode
 */
typedef enum {ROOT=0, MIN_MAX=1, INFLECTION=2} mode;

/*
 * Signature for solver model functions
 */
typedef void (*model)(series, series);

/*
 * Initialize file scoped temporary storage
 */
void ad_tempvars (void);

/*
 * Selects a jet for automatic diffentiation, or not
 */
void set_ad_status (series jet, ad_status s);

/*
 * Prints a Taylor coefficient jet to order n
 */
void jet_output (series jet, long n, char* f_colour, char *fk_colour);

/*
 * Applies factorials to convert Taylor coefficients to actual derivative values
 */
void jet_to_derivs (series jet, long n);

/*
 * Prints a derivative jet to order n
 */
void derivative_output (series jet, long n, char* f_colour, char *fk_colour);

/*
 * Finds a root of fn(f, x) by Newton's method, where f and x are Taylor Series
 */
void ad_newton (model m, series f, series x, int max_it, mpfr_t f_tol, mpfr_t x_tol, mode degree);

series ad_set (series B, series A);

/*
 * Scales Taylor Series U by a factor a, result stored in jet S
 */
series ad_scale (series S, series U, mpfr_t a);

/*
 * Calculates Taylor Series for the sum of U and V, result stored in jet P
 */
series ad_plus (series P, series V, series U);

/*
 * Calculates Taylor Series for the difference of U and V, result stored in jet M
 */
series ad_minus (series M, series V, series U);

/*
 * Calculates Taylor Series for the negative of U, result stored in jet M
 */
series ad_neg (series M, series U);

/*
 * Calculates Taylor Series for the absolute value of U, result stored in jet A
 */
series ad_abs (series A, series U);

/*
 * Calculates Taylor Series for U * V, result stored in jet P
 */
series ad_prod (series P, series V, series U);

/*
 * Calculates Taylor Series for U / V, result stored in jet Q
 */
series ad_quot (series Q, series U, series V);

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
tuple ad_sin_cos (series S, series C, series U, geometry g);

/*
 * Calculates Taylor Series for the trigonometric tangent and squared secant of U, results stored in jets T and S2
 */
tuple ad_tan_sec2 (series T, series S2, series U, geometry g);

/*
 * Calculates Taylor Series for U^a, results stored in jet P
 */
series ad_pwr (series P, series U, double a);

/*
 * Calculates Taylor Series for the natural logarithm of U, results stored in jet L
 */
series ad_ln (series L, series U);
