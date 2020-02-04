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
    mpfr_t *a;
    mpfr_t *b;
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
 * Returns a pointer to kth element of the absolute value of U, result stored and returned in variable A, NO JET STORAGE
 */
mpfr_t *t_abs (mpfr_t *a, mpfr_t *u, int k);

/*
 * Cauchy product for C = A.B
 *
 *  let c(x) = sum{k=0->inf} C(k) (x - a)^k
 *
 *   if c(x) = a(x) b(x)
 *
 * then c(x) = (sum{k=0->inf} A(k) (x - a)^k) (sum{k=0->inf} B(k) (x - a)^k)
 *
 *           = sum{k=0->inf} ( sum{j=0->k} A[j]B[k - j] ) (x - a)^k
 *
 *  ==> C[k] = sum{j=0->k} A[j]B[k - j]
 */

/*
 * Returns a pointer to kth element of the square of U, result stored and returned in variable S, NO JET STORAGE
 *
 *  S = U.U
 *
 *  S = sum{j=0->k} U[j].U[k-j]
 *
 *  or, from symmetry:
 *
 *  S = sum{j=0->(k-1)/2} U[j].U[k-j]             if k odd
 *  S = sum{j=0->(k-2)/2} U[j].U[k-j] + U[k/2]^2  if k even
 */
mpfr_t *t_sqr (mpfr_t *S, mpfr_t *U, int k);

/*
 * Returns a pointer to kth element of the product of U and V, result stored in variable P, NO JET STORAGE
 *
 *  P = U.V
 *
 *  P = sum{j=0->k} U[j].V[k-j]
 */
mpfr_t *t_prod (mpfr_t *P, mpfr_t *U, mpfr_t *V, int k);

/*
 * Returns a pointer to kth element of U / V, results accumulated in jet Q, DOMAIN RESTRICTION v[0] != 0.0
 *
 *    Q = U / V ==> U = Q.V
 *
 *                U[k] = sum{j=0->k} Q[j].V[k-j]
 *
 *                     = sum{j=0->k-1} Q[j].V[k-j] + Q[k].V[0]
 *
 *                Q[k] = (U[k] - sum{j=0->k-1} Q[j].V[k-j]) / V[0]
 */
mpfr_t *t_quot (mpfr_t *Q, mpfr_t *U, mpfr_t *V, int k);

/*
 * Returns a pointer to kth element of the square root of U, results accumulated in jet R, DOMAIN RESTRICTION U[0] > 0.0
 *
 *    U = R.R
 *
 * U[k] = sum{j=0->k} R[j].R[k - j]
 *
 *      = sum{j=1->k-1} R[j].R[k - j] + 2.R[k].R[0]
 *
 * R[k] = (U[k] - sum{j=1->k-1} R[j].R[k-j]) / (2.R[0])
 *
 * or, from symmetry:
 *
 * R[k] = (U[k] - sum{j=1->(k-1)/2} R[j].R[k-j]) / (2.0.R[0])             if k odd
 * R[k] = (U[k] - sum{j=1->(k-2)/2} R[j].R[k-j] - R[k/2]^2) / (2.0.R[0])  if k even
 */
mpfr_t *t_sqrt (mpfr_t *R, mpfr_t *U, int k);

/*
 * Applying the chain rule for the derivative of a composed function f(u) creates another Cauchy product:
 *
 *          F' = (df/du).U'
 *             =       H.U'
 *
 * Using F'[k] = (k+1)F[k+1]  ==>  F'[k-1] = kF[k], we can replace F' with F, and U' with U as follows:
 *
 * from product rule above, (note that F' and U' have one fewer elements than F and U)
 *
 *     F'[k-1] = sum{j=0->k-1} H[j].U'[k-1-j]
 *
 *       kF[k] = sum{j=0->k-1} H[j].(k-j)U[k-j]
 *
 *        F[k] = sum{j=0->k-1} H[j].(k-j)U[k-j]/k
 */

/*
 * Returns a pointer to kth element of the exponential of U, results accumulated in jet E
 *
 *      E' = E.U'
 *
 *    E[k] = sum{j=0->k-1} E[j].(k-j)U[k-j]/k
 */
mpfr_t *t_exp (mpfr_t *E, mpfr_t *U, int k, mpfr_t *tmp);

/*
 * Returns a pair of pointers to kth elements of the sine and cosine of U, results accumulated in jets S and C
 *
 *      S' =       C.U'
 *      C' = (+/-) S.U'   (+ for cosh(g=HYP), - for cos(g=TRIG))
 *
 *    S[k] = sum{j=0->k-1}       C[j].(k-j)U[k-j]/k
 *    C[k] = sum{j=0->k-1} (+/-) S[j].(k-j)U[k-j]/k
 */
tuple t_sin_cos (mpfr_t *S, mpfr_t *C, mpfr_t *U, int k, mpfr_t *tmp, geometry g);

/*
 * Returns a pair of pointers to kth elements of the tangent and squared secant of U, results accumulated in jets T and S2
 *
 *      T' =       S2.U'
 *     S2' = (+/-)2 T.T'   (+ for sec^2(g=TRIG), - for sech^2(g=HYP))
 *
 *    T[k] = sum{j=0->k-1}       S2[j].(k-j)U[k-j]/k
 *   S2[k] = sum{j=0->k-1} (+/-)2 T[j].(k-j)T[k-j]/k
 */
tuple t_tan_sec2 (mpfr_t *T, mpfr_t *S2, mpfr_t *U, int k, mpfr_t *tmp, geometry g);

/*
 * Returns a pointer to kth element of P = U^a (where a is scalar), results accumulated in jet P, DOMAIN RESTRICTION U[0] > 0.0
 *
 *                    P'= U^a' = a U^(a-1).U'
 *                      U.U^a' = a U^a.U'
 *                        U.P' = a P.U'
 *
 *    sum{j=0->k} U[k-j].jP[j] =  a sum{j=0->k-1} P[j].(k-j)U[k-j]
 *
 *                  U[0].kP[k] =  a sum{j=0->k-1} P[j].(k-j)U[k-j]   - sum{j=0->k-1} U[k-j].jP[j]
 *
 *                        P[k] = (a sum{j=0->k-1} P[j].(k-j)U[k-j]/k - sum{j=0->k-1} U[k-j].jP[j]/k) / U[0]
 */
mpfr_t *t_pwr (mpfr_t *P, mpfr_t *U, double a, int k, mpfr_t *tmpA, mpfr_t *tmpB);

/*
 * Returns a pointer to kth element of the natural logarithm of U, result accumulated in jet L, DOMAIN RESTRICTION U[0] > 0.0
 *
 *     L' = U' / U ==> U' = U.L'
 *
 *                 k.U[k] = sum{j=0->k-1} U[j].(k-j)L[k-j]
 *
 *                        = sum{j=1->k-1} U[j].(k-j)L[k-j] + U[0].kL[k]
 *
 *                   L[k] = (U[k] - sum{j=1->k-1} U[j].(k-j)L[k-j]/k) / U[0]
 */
mpfr_t *t_ln (mpfr_t *L, mpfr_t *U, int k, mpfr_t *tmp);

