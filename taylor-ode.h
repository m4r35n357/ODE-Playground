/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

/*
 * Diffentiate real numbers from jets
 */
typedef long double real;

typedef long double *series;

/*
 * Prints an index column, and x, y, z columns, into a single line
 */
void t_output (long dp, real x, real y, real z, real t, char *x_label, char *y_label, char *z_label) ;

/*
 * Sets output precision, integrator parameters and initial conditions from the command line arguments (1 to 4)
 */
void t_control (char **argv, long *dp, long *n, real *h, long *nsteps, real *x, real *y, real *z);

/*
 * Bulk set ODE parameters from the command line arguments (5 onwards)
 */
void t_params (char **argv, int count, ...);

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
 * Returns kth element of the exponential of U, results accumulated in jet E
 */
real t_exp (series E, series U, int k);

/*
 * Selects either a trigonometric or hyperbolic version of the function
 */
typedef enum {TRIG, HYP} geometry;

/*
 * For returning combined recurrence values
 */
typedef struct {
    real a;
    real b;
} pair;

/*
 * Returns kth elements of both sine and cosine of U, results accumulated in jets S and C
 */
pair t_sin_cos (series S, series C, series U, int k, geometry g);

/*
 * Returns kth elements of both tangent and squared secant of U, results accumulated in jets T and S2
 */
pair t_tan_sec2 (series T, series S2, series U, int k, geometry g);

/*
 * Function signatures for ODE parameters and intermediate variables, defined in client code
 */
typedef void *(*tsm_params)(int, char **, long);

typedef void *(*tsm_inters)(long);

/*
 * For returning x, y, z values
 */
typedef struct {
    real x;
    real y;
    real z;
} components;

/*
 * 3D (x, y, z) ODE function signatures, defined in client code
 */
typedef components (*tsm_model)(series, series, series, void *, void *, int);

/*
 * Integrator signatures
 */
void tsm (int argc, char **argv, tsm_model ode, tsm_params get_p, tsm_inters get_i);
