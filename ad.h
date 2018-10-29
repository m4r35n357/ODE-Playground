
/*
 * Interface for performing high order Automatic Differentiation
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

/*
 * Signature for solver model functions
 */
typedef void (*model)(mpfr_t *, mpfr_t *, int);

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
 * Finds a root of fn(f, x) by bisection, where f and x are Taylor Series
 */
int ad_bisect (model m, mpfr_t *xa, mpfr_t *xb, int max_it, mpfr_t f_tol, mpfr_t x_tol, mpfr_t *xc, mpfr_t *fa, mpfr_t *fc);

/*
 * Finds a root of fn(f, x) by Newton's method, where f and x are Taylor Series
 */
int ad_newton (model m, mpfr_t *f, mpfr_t *x, int max_it, mpfr_t f_tol, mpfr_t x_tol);

/*
 * Finds a root of fn(f, x) by Householder's method, where f and x are Taylor Series
 */
int ad_householder (model m, mpfr_t *f, mpfr_t *x, long n, int max_it, mpfr_t f_tol, mpfr_t x_tol, mpfr_t *f_reciprocal, mpfr_t *w1);

/*
 * Scales Taylor Series U by a factor a, result stored in jet S
 */
void ad_scale (mpfr_t *S, mpfr_t *U, mpfr_t a, int n);

/*
 * Calculates Taylor Series for the sum of U and V, result stored in jet P
 */
void ad_plus (mpfr_t *P, mpfr_t *V, mpfr_t *U, int n);

/*
 * Calculates Taylor Series for the difference of U and V, result stored in jet P
 */
void ad_minus (mpfr_t *P, mpfr_t *V, mpfr_t *U, int n);

/*
 * Calculates Taylor Series for the square of U, result stored in jet S
 */
void ad_square (mpfr_t *S, mpfr_t *U, int n);

/*
 * Calculates Taylor Series for U * V, result stored in jet P
 */
void ad_product(mpfr_t *P, mpfr_t *V, mpfr_t *U, int n);

/*
 * Calculates Taylor Series for U / V, result stored in jet Q
 */
void ad_quotient (mpfr_t *Q, mpfr_t *U, mpfr_t *V, int n);

/*
 * Calculates Taylor Series for the exponential of U, results stored in jet E
 */
void ad_exp (mpfr_t *E, mpfr_t *U, int n, mpfr_t *tmp);

/*
 * Calculates Taylor Series for the trigonometric sine and cosine of U, results stored in jets S and C
 */
void ad_sin_cos (mpfr_t *S, mpfr_t *C, mpfr_t *U, int n, mpfr_t *tmp);

/*
 * Calculates Taylor Series for the hyperbolic sine and cosine of U, results stored in jets S and C
 */
void ad_sinh_cosh (mpfr_t *S, mpfr_t *C, mpfr_t *U, int n, mpfr_t *tmp);

/*
 * Calculates Taylor Series for the trigonometric tangent and squared secant of U, results stored in jets T and S2
 */
void ad_tan_sec2 (mpfr_t *T, mpfr_t *S2, mpfr_t *U, int n, mpfr_t *tmp);

/*
 * Calculates Taylor Series for the hyperbolic tangent and squared secant of U, results stored in jets T and S2
 */
void ad_tanh_sech2 (mpfr_t *T, mpfr_t *S2, mpfr_t *U, int n, mpfr_t *tmp);

/*
 * Calculates Taylor Series for U^a, results stored in jet P
 */
void ad_power (mpfr_t *P, mpfr_t *U, mpfr_t a, int n, mpfr_t *tmp1, mpfr_t *tmp2);

/*
 * Calculates Taylor Series for the natural logarithm of U, results stored in jet L
 */
void ad_ln (mpfr_t *L, mpfr_t *U, int n, mpfr_t *tmp);
