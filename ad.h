
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
typedef enum {ROOT___=0, MIN_MAX=1, INFLECT=2} mode;

/*
 * Signature for solver model functions
 */
typedef void (*model)(mpfr_t *, mpfr_t *, int);

/*
 * Selects a jet for automatic diffentiation, or not
 */
void set_ad_status (mpfr_t *jet, ad_status s);

/*
 * Prints a Taylor coefficient jet to order n
 */
void jet_output (mpfr_t *jet, long n, char* f_colour, char *fk_colour);

/*
 * Applies factorials to convert Taylor coefficients to actual derivative values
 */
void jet_to_derivs (mpfr_t *jet, long n);

/*
 * Prints a derivative jet to order n
 */
void derivative_output (mpfr_t *jet, long n, char* f_colour, char *fk_colour);

/*
 * Finds a root of fn(f, x) by Newton's method, where f and x are Taylor Series
 */
void ad_newton (model m, mpfr_t *f, mpfr_t *x, int max_it, mpfr_t f_tol, mpfr_t x_tol, mode degree, mpfr_t *delta);

mpfr_t *ad_set (mpfr_t *B, mpfr_t *A, int n);

/*
 * Scales Taylor Series U by a factor a, result stored in jet S
 */
mpfr_t *ad_scale (mpfr_t *S, mpfr_t *U, mpfr_t a, int n);

/*
 * Calculates Taylor Series for the sum of U and V, result stored in jet P
 */
mpfr_t *ad_plus (mpfr_t *P, mpfr_t *V, mpfr_t *U, int n);

/*
 * Calculates Taylor Series for the difference of U and V, result stored in jet P
 */
mpfr_t *ad_minus (mpfr_t *P, mpfr_t *V, mpfr_t *U, int n);

/*
 * Calculates Taylor Series for the absolute value of U, result stored in jet A
 */
mpfr_t *ad_abs (mpfr_t *A, mpfr_t *U, int n);

/*
 * Calculates Taylor Series for U * V, result stored in jet P
 */
mpfr_t *ad_prod (mpfr_t *P, mpfr_t *V, mpfr_t *U, int n);

/*
 * Calculates Taylor Series for U / V, result stored in jet Q
 */
mpfr_t *ad_quot (mpfr_t *Q, mpfr_t *U, mpfr_t *V, int n);

/*
 * Calculates Taylor Series for the square of U, result stored in jet S
 */
mpfr_t *ad_sqr (mpfr_t *S, mpfr_t *U, int n);

/*
 * Calculates Taylor Series for the square root of U, result stored in jet R
 */
mpfr_t *ad_sqrt (mpfr_t *R, mpfr_t *U, int n);

/*
 * Calculates Taylor Series for the exponential of U, results stored in jet E
 */
mpfr_t *ad_exp (mpfr_t *E, mpfr_t *U, int n);

/*
 * Calculates Taylor Series for the trigonometric sine and cosine of U, results stored in jets S and C
 */
tuple ad_sin_cos (mpfr_t *S, mpfr_t *C, mpfr_t *U, int n, geometry g);

/*
 * Calculates Taylor Series for the trigonometric tangent and squared secant of U, results stored in jets T and S2
 */
tuple ad_tan_sec2 (mpfr_t *T, mpfr_t *S2, mpfr_t *U, int n, geometry g);

/*
 * Calculates Taylor Series for U^a, results stored in jet P
 */
mpfr_t *ad_pwr (mpfr_t *P, mpfr_t *U, double a, int n);

/*
 * Calculates Taylor Series for the natural logarithm of U, results stored in jet L
 */
mpfr_t *ad_ln (mpfr_t *L, mpfr_t *U, int n);
