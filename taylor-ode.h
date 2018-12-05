
/*
 * Interface for solving systems of Ordinary Differentrial Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

/*
 * The numerical base for string IO conversions
 */
const int BASE;

/*
 * Global rounding strategy for MPFR
 */
const mpfr_rnd_t RND;

/*
 * Selects either a trigonometric or hyperbolic version of the function
 */
typedef enum {TRIG, HYP} geometry;

/*
 * Prints an index column, and any number of data columns, into a single line
 */
void t_line_output (mpfr_t t, int item_count, ...);

/*
 * Sets the order, step size, and number of steps for the integration
 */
void t_stepper (char **argv, long *n, mpfr_t *t, mpfr_t *h, long *nsteps);

/*
 * Creates an initialized derivative jet of the specified size, with no values set
 */
mpfr_t *t_jet (int size);

/*
 * Creates a derivative jet of the specified size, with element zero set to value and the rest zeroed
 */
mpfr_t *t_jet_c (int size, mpfr_t value);

/*
 * Sums a Taylor series safely and efficiently
 */
void t_horner (mpfr_t *sum, const mpfr_t *jet, long n, mpfr_t h);

/*
 * Calculates kth element of the square of U, result stored in variable S
 */
void t_square (mpfr_t *S, const mpfr_t *U, int k);

/*
 * Calculates kth element of the product of U and V, result stored in variable P
 */
void t_product (mpfr_t *P, const mpfr_t *U, const mpfr_t *V, int k);

/*
 * Calculates kth element of U / V, results stored in jet Q
 */
void t_quotient (mpfr_t *Q, const mpfr_t *U, const mpfr_t *V, int k);

/*
 * Calculates kth element of the square root of U, results stored in jet R
 */
void t_sqrt (mpfr_t *R, const mpfr_t *U, int k);

/*
 * Calculates kth element of the exponential of U, results stored in jet E
 */
void t_exp (mpfr_t *E, const mpfr_t *U, int k, mpfr_t *tmp);

/*
 * Calculates kth elements of the sine and cosine of U, results stored in jets S and C
 */
void t_sin_cos(mpfr_t *S, mpfr_t *C, const mpfr_t *U, int k, mpfr_t *tmp, geometry g);

/*
 * Calculates Taylor Series for the tangent and squared secant of U, results stored in jets T and S2
 */
void t_tan_sec2 (mpfr_t *T, mpfr_t *S2, const mpfr_t *U, int k, mpfr_t *tmp, geometry g);

/*
 * Calculates kth element of U^a, results stored in jet P
 */
void t_power (mpfr_t *P, const mpfr_t *U, mpfr_t a, int k, mpfr_t *tmp1, mpfr_t *tmp2);

/*
 * Calculates Taylor Series for U * V, result stored in jet P
 */
void t_ln (mpfr_t *L, const mpfr_t *U, int k, mpfr_t *tmp);

