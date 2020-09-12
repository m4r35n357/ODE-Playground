/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

/*
 * Isolate real numbers from jets
 */
typedef long double real;
typedef long double *series;

/*
 * The numerical base for string IO conversions
 */
const int BASE;

/*
 * For returning combined recurrence values
 */
typedef struct {
    real a;
    real b;
} pair;

/*
 * Selects either a trigonometric or hyperbolic version of the function
 */
typedef enum {TRIG, HYP} geometry;

/*
 * Prints an index column, and x, y, z columns, into a single line
 */
void t_xyz_output (real x, real y, real z, real t);

/*
 * Sets the order, step size, and real of steps for the integration from the command line arguments (1 to 4)
 */
void t_stepper (char **argv, long *n, real *h, long *nsteps);

/*
 * Bulk set initial conditions and ODE parameters from the command line arguments (5 onwards)
 */
void t_args (char **argv, int count, ...);

/*
 * Creates an initialized jet of the specified size, with no values set
 */
series t_jet (int size);

/*
 * Returns a jet of the specified size, with element zero set to value and the rest zeroed (represents a constant in an ODE)
 */
series t_jet_c (int size, real value);

/*
 * Sums a Taylor series safely and efficiently.  Result is BOTH returned and stored in jet[0]
 */
real t_horner (series jet, int n, real h);

/*
 * Returns kth element of the absolute value of U, NO JET STORAGE
 */
real t_abs (series u, int k);

/*
 * Returns kth element of the square of U, NO JET STORAGE
 */
real t_sqr (series U, int k);

/*
 * Returns kth element of the product of U and V, NO JET STORAGE
 */
real t_prod (series U, series V, int k);

/*
 * Returns kth element of U / V, results accumulated in jet Q, DOMAIN RESTRICTION v[0] != 0.0
 */
real t_quot (series Q, series U, series V, int k);

/*
 * Returns kth element of 1 / V, results accumulated in jet I, DOMAIN RESTRICTION v[0] != 0.0
 */
real t_inv (series I, series V, int k);

/*
 * Returns kth element of the square root of U, results accumulated in jet R, DOMAIN RESTRICTION U[0] > 0.0
 */
real t_sqrt (series R, series U, int k);

/*
 * Returns kth element of the exponential of U, results accumulated in jet E
 */
real t_exp (series E, series U, int k);

/*
 * Returns struct containing kth elements of both sine and cosine of U, results accumulated in jets S and C
 */
pair t_sin_cos (series S, series C, series U, int k, geometry g);

/*
 * Returns struct containing kth elements of both tangent and squared secant of U, results accumulated in jets T and S2
 */
pair t_tan_sec2 (series T, series S2, series U, int k, geometry g);

/*
 * Returns kth element of P = U^a (where a is scalar), results accumulated in jet P, DOMAIN RESTRICTION U[0] > 0.0
 */
real t_pwr (series P, series U, real a, int k);

/*
 * Returns kth element of the natural logarithm of U, result accumulated in jet L, DOMAIN RESTRICTION U[0] > 0.0
 */
real t_ln (series L, series U, int k);

