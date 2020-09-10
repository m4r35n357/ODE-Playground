/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

/*
 * The numerical base for string IO conversions
 */
const int BASE;

/*
 * For returning combined recurrence values
 */
typedef struct {
    long double *a;
    long double *b;
} tuple;

/*
 * Selects either a trigonometric or hyperbolic version of the function
 */
typedef enum {TRIG, HYP} geometry;

/*
 * Prints an index column, and x, y, z columns, into a single line
 */
void t_xyz_output (long double x, long double y, long double z, long double t);

/*
 * Sets the order, step size, and number of steps for the integration from the command line arguments (1 to 4)
 */
void t_stepper (char **argv, long *n, long double *h, long *nsteps);

/*
 * Bulk set initial conditions and ODE parameters from the command line arguments (5 onwards)
 */
void t_args (char **argv, int count, ...);

/*
 * Creates an initialized jet of the specified size, with no values set
 */
long double *t_jet (int size);

/*
 * Returns a jet of the specified size, with element zero set to value and the rest zeroed (represents a constant in an ODE)
 */
long double *t_jet_c (int size, long double value);

/*
 * Sums a Taylor series safely and efficiently.  Result is BOTH returned and stored in jet[0]
 */
long double t_horner (long double *jet, int n, long double h);

/*
 * Returns kth element of the absolute value of U, NO JET STORAGE
 */
long double t_abs (long double *u, int k);

/*
 * Returns kth element of the square of U, NO JET STORAGE
 */
long double t_sqr (long double *U, int k);

/*
 * Returns kth element of the product of U and V, NO JET STORAGE
 */
long double t_prod (long double *U, long double *V, int k);

/*
 * Returns kth element of U / V, results accumulated in jet Q, DOMAIN RESTRICTION v[0] != 0.0
 */
long double t_quot (long double *Q, long double *U, long double *V, int k);

/*
 * Returns kth element of 1 / V, results accumulated in jet I, DOMAIN RESTRICTION v[0] != 0.0
 */
long double t_inv (long double *I, long double *V, int k);

/*
 * Returns kth element of the square root of U, results accumulated in jet R, DOMAIN RESTRICTION U[0] > 0.0
 */
long double t_sqrt (long double *R, long double *U, int k);

/*
 * Returns kth element of the exponential of U, results accumulated in jet E
 */
long double t_exp (long double *E, long double *U, int k);

/*
 * Returns struct of pointers to kth elements of both sine and cosine of U, results accumulated in jets S and C
 */
tuple t_sin_cos (long double *S, long double *C, long double *U, int k, geometry g);

/*
 * Returns struct of pointers to kth elements of both tangent and squared secant of U, results accumulated in jets T and S2
 */
tuple t_tan_sec2 (long double *T, long double *S2, long double *U, int k, geometry g);

/*
 * Returns kth element of P = U^a (where a is scalar), results accumulated in jet P, DOMAIN RESTRICTION U[0] > 0.0
 */
long double t_pwr (long double *P, long double *U, long double a, int k);

/*
 * Returns kth element of the natural logarithm of U, result accumulated in jet L, DOMAIN RESTRICTION U[0] > 0.0
 */
long double t_ln (long double *L, long double *U, int k);

