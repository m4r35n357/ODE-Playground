/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018,2019 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
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
 * For returning "paired" recurrence values
 */
struct Tuple {
    mpfr_t *const a;
    mpfr_t *const b;
};

/*
 * Selects either a trigonometric or hyperbolic version of the function
 */
typedef enum {TRIG, HYP} geometry;

/*
 * Set an MPFR variable from a numbered string argument
 */
void t_arg (char **argv, int arg_no, mpfr_t *variable);

/*
 * Prints an index column, and x, y, z columns, into a single line
 */
void t_xyz_output (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t t);

/*
 * Sets the order, step size, and number of steps for the integration
 */
void t_stepper (char **argv, long *n, mpfr_t *t, mpfr_t *h, long *nsteps);

/*
 * Creates an initialized jet of the specified size, with no values set
 */
mpfr_t *t_jet (int size);

/*
 * Returns a jet of the specified size, with element zero set to value and the rest zeroed
 */
mpfr_t *t_jet_c (int size, mpfr_t value);

/*
 * Sums a Taylor series safely and efficiently
 */
void t_horner (mpfr_t *sum, mpfr_t *jet, int n, mpfr_t h);

/*
 * Returns kth element of the absolute value of U, result stored and returned in variable A
 */
mpfr_t *t_abs (mpfr_t *a, mpfr_t *u, int k);

/*
 * Returns kth element of the square of U, result stored and returned in variable S
 */
mpfr_t *t_sqr (mpfr_t *S, mpfr_t *U, int k);

/*
 * Returns kth element of the product of U and V, result stored in variable P
 */
mpfr_t *t_prod (mpfr_t *P, mpfr_t *U, mpfr_t *V, int k);

/*
 * Returns kth element of U / V, results stored in jet Q
 */
mpfr_t *t_quot (mpfr_t *Q, mpfr_t *U, mpfr_t *V, int k);

/*
 * Returns kth element of the square root of U, results stored in jet R
 */
mpfr_t *t_sqrt (mpfr_t *R, mpfr_t *U, int k);

/*
 * Returns kth element of U^a, results stored in jet P
 */
mpfr_t *t_pwr (mpfr_t *P, mpfr_t *U, mpfr_t a, int k, mpfr_t *tmp);

/*
 * Returns kth element of the exponential of U, results stored in jet E
 */
mpfr_t *t_exp (mpfr_t *E, mpfr_t *U, int k, mpfr_t *tmp);

/*
 * Returns kth element of the natural logarithm of U, result stored in jet L
 */
mpfr_t *t_ln (mpfr_t *L, mpfr_t *U, int k, mpfr_t *tmp);

/*
 * Returns a struct containing kth elements of the sine and cosine of U, results stored in jets S and C
 */
struct Tuple t_sin_cos(mpfr_t *S, mpfr_t *C, mpfr_t *U, int k, mpfr_t *tmp, geometry g);

/*
 * Returns a struct containing kth elements of the tangent and squared secant of U, results stored in jets T and S2
 */
struct Tuple t_tan_sec2 (mpfr_t *T, mpfr_t *S2, mpfr_t *U, int k, mpfr_t *tmp, geometry g);

