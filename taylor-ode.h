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
 * Global rounding strategy for MPFR
 */
const mpfr_rnd_t RND;

/*
 * For returning "paired" recurrence values
 */
typedef struct {
    mpfr_t *const a;
    mpfr_t *const b;
} tuple;

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
 * Returns a jet of the specified size, with element zero set to value and the rest zeroed (represents a constant in an ODE)
 */
mpfr_t *t_jet_c (int size, mpfr_t value);

/*
 * Sums a Taylor series safely and efficiently
 */
void t_horner (mpfr_t *sum, mpfr_t *jet, int n, mpfr_t h);

/*
 * Returns kth element of the absolute value of U, result stored and returned in variable A, NO JET STORAGE
 */
mpfr_t *t_abs (mpfr_t *a, mpfr_t *u, int k);

/*
 * Product rule for c = a.b
 *
 *  c[k] = sum{j = 0 to k} a[j].b[k - j]
 */

/*
 * Returns kth element of the square of U, result stored and returned in variable S, NO JET STORAGE
 *
 *  s = u.u
 *
 * s = sum{j = 0 to (k - 1) / 2} u[j].u[k - j]               if k odd
 * s = sum{j = 0 to (k - 2) / 2} u[j].u[k - j] + u[k / 2]^2  if k even
 */
mpfr_t *t_sqr (mpfr_t *S, mpfr_t *U, int k);

/*
 * Returns kth element of the product of U and V, result stored in variable P, NO JET STORAGE
 *
 *  p = u.v
 *
 * p = sum{j = 0 to k} u[j].v[k - j]
 */
mpfr_t *t_prod (mpfr_t *P, mpfr_t *U, mpfr_t *V, int k);

/*
 * Returns kth element of U / V, results stored in jet Q, DOMAIN RESTRICTION v[0] != 0.0
 *
 *  u = q.v
 *
 * q[k] = (u[k] - sum{j = 0 to k - 1} q[j].v[k - j]) / v[0]
 */
mpfr_t *t_quot (mpfr_t *Q, mpfr_t *U, mpfr_t *V, int k);

/*
 * Returns kth element of the square root of U, results stored in jet R, DOMAIN RESTRICTION U[0] > 0.0
 *
 *  u = r.r
 *
 * r[k] = (u[k] - sum{j = 0 to (k - 1) / 2} r[j].r[k - j]) / (2.0 * u[0])               if k odd
 * r[k] = (u[k] - sum{j = 0 to (k - 2) / 2} r[j].r[k - j] - r[k / 2]^2) / (2.0 * u[0])  if k even
 */
mpfr_t *t_sqrt (mpfr_t *R, mpfr_t *U, int k);

/*
 * Chain rule for derivative of composed function f(u) creates another product rule:
 *
 *  f'(u) = (df/du).u'
 *
 * From product rule above, (note that f' and u' have one fewer elements than f and u)
 *
 *  f'[k - 1] = sum{j = 0 to k - 1} (df/du)[j].(k - 1 - j).u'[k - 1 - j]
 *
 * Using f'[k] = (k + 1).f[k + 1], or f'[k - 1] = k.f[k], we can replace f' with f, and u' with u
 *
 *  -> k.f[k] = sum{j = 0 to k - 1} (df/du)[j].(k - j).u[k - j]
 */

/*
 * Returns kth element of the exponential of U, results stored in jet E
 *
 * exp'(u) = exp(u).u'
 *
 * exp(u)[k] = sum{j = 0 to k - 1} exp(u)[j].(k - j).u[k - j] / k
 */
mpfr_t *t_exp (mpfr_t *E, mpfr_t *U, int k, mpfr_t *tmp);

/*
 * Returns a struct containing kth elements of the sine and cosine of U, results stored in jets S and C
 *
 *  sinh'(u) = cosh(u).u', sin'(u) =  cos(u).u'
 *  cosh'(u) = sinh(u).u', cos'(u) = -sin(u).u'
 *
 *  sin(u)[k] = sum{j = 0 to k - 1}  cos(u)[j].(k - j).u[k - j] / k
 *  cos(u)[k] = sum{j = 0 to k - 1} -sin(u)[j].(k - j).u[k - j] / k
 *
 * sinh(u)[k] = sum{j = 0 to k - 1} cosh(u)[j].(k - j).u[k - j] / k
 * cosh(u)[k] = sum{j = 0 to k - 1} sinh(u)[j].(k - j).u[k - j] / k
 */
tuple t_sin_cos (mpfr_t *S, mpfr_t *C, mpfr_t *U, int k, mpfr_t *tmp, geometry g);

/*
 * Returns a struct containing kth elements of the tangent and squared secant of U, results stored in jets T and S2
 *
 *    tanh'(u) =  sech^2(u).u',         tan'(u) = sec^2(u).u'
 *  sech^2'(u) = -2.tanh(u).tanh'(u), sec^2'(u) = 2.tan(u).tan'(u)
 *
 *    tan(u)[k] = sum{j = 0 to k - 1} sec^2(u)[j].(k - j).u[k - j] / k
 *  sec^2(u)[k] = sum{j = 0 to k - 1} 2.tan(u)[j].(k - j).tan(u)[k - j] / k
 *
 *   tanh(u)[k] = sum{j = 0 to k - 1}  sech^2(u)[j].(k - j).u[k - j] / k
 * sech^2(u)[k] = sum{j = 0 to k - 1} -2.tanh(u)[j].(k - j).tanh(u)[k - j] / k
 */
tuple t_tan_sec2 (mpfr_t *T, mpfr_t *S2, mpfr_t *U, int k, mpfr_t *tmp, geometry g);

/*
 * Returns kth element of U^a (where a is scalar), results stored in jet P, DOMAIN RESTRICTION U[0] > 0.0
 *
 *  u^a'.u = a.u^a.u'
 *       0 = a.u^a.u' - u^a'.u
 *
 *           0 = sum{j = 0 to k - 1} a.p[j].(k - j).u[k - j] - sum{j = 0 to k} j.p[j].u[k - j]
 * k.p[k].u[0] = sum{j = 0 to k - 1} a.p[j].(k - j).u[k - j] - sum{j = 0 to k - 1} j.p[j].u[k - j]
 *
 *        p[k] = sum{j = 0 to k - 1} (a.(k - j) - j).p[j].u[k - j] / (k * u[0])
 */
mpfr_t *t_pwr (mpfr_t *P, mpfr_t *U, mpfr_t a, int k, mpfr_t *tmp);

/*
 * Returns kth element of the natural logarithm of U, result stored in jet L, DOMAIN RESTRICTION U[0] > 0.0
 *
 *  u' = u.ln'(u)
 *
 * k.u[k] = sum{j = 0 to k - 1} u[j].(k - j).ln[k - j]
 *        = sum{j = 1 to k - 1} u[j].(k - j).ln[k - j] + u[0].k.ln[k]
 *
 * ln[k] = (u[k] - sum{j = 1 to k - 1} u[j] * (k - j) * ln[k - j] / k) / u[0]
 */
mpfr_t *t_ln (mpfr_t *L, mpfr_t *U, int k, mpfr_t *tmp);

