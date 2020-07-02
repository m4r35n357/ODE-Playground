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
 * Wrap coefficient "jets" in a structure
 */
typedef struct {
    mpfr_t *jet;
    int size;
} series;

/*
 * For returning "paired" recurrence values
 */
typedef struct {
    mpfr_t *a;
    mpfr_t *b;
    int size;
} tuple;

/*
 * Sign of a number
 */
typedef enum {POS, NEG} sign;

/*
 * Selects the trigonometric or hyperbolic version of a function
 */
typedef enum {TRIG, HYP} geometry;

/*
 * Pre-allocate some MPFR variaboles
 */
void t_tempvars (void);

/*
 * Prints x, y, z, t values in a single line of output
 */
void t_output (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t h, long step);

/*
 * Sets the (precision,) order, step size, and number of steps for the integration from the command line arguments (1 to 5)
 */
void t_stepper (char **argv, long *n, mpfr_t *h, long *nsteps);

/*
 * Bulk set initial coordinate values and ODE parameters from the command line arguments (6 onwards)
 */
void t_args (char **argv, int count, ...);

/*
 * Creates a zeroed Taylor Series of the specified size
 */
series t_series (int size);

/*
 * The Taylor Series Method (TSM) in brief:
 *
 *              x(t0 + h) = x(t) = sum{k=0->inf} X[k] h^k,    where X[0] = x(t0), X[k] = (d/dt)^k x(t0) / k! and h = t - t0
 *
 *                         x'(t) = sum{k=0->inf} X'[k] h^k,   where x'(t) = dx/dt, the ODE equations   (A)
 *
 *                     d/dt x(t) = sum{k=1->inf} k X[k] h^(k-1)
 *
 *                               = sum{k=0->inf} (k+1) X[k+1] h^k                                      (B)
 *
 * Comparing (A) and (B),  X'[k] = (k+1) X[k+1]    *** this is THE IDENTITY (also used in recurrences below) ***
 *
 *                   ==>  X[k+1] = X'[k] / (k + 1)
 *
 * 1. Build up a coefficient jet for each x(t) using the ODEs, x'(t), and Taylor recurrences where needed, then divide by k + 1
 *
 * 2. Apply Horner's method to each jet to calculate the next set of coordinates
 */

/*
 * Calculate next coefficient in jet
 */
void t_next (series S, mpfr_t dot, int k, sign sgn);

/*
 * Evaluate a Taylor series safely and efficiently
 */
mpfr_t *t_horner (series S, mpfr_t h);

/*
 * Returns a pointer to kth element of the absolute value of U, result stored and returned in variable A, NO JET STORAGE
 */
mpfr_t *t_abs (series U, int k);

/*
 * Cauchy product for C = A.B
 *
 *   if c(t) = a(t) b(t)
 *
 * then c(t) = sum{k=0->inf} C(k) h^k
 *
 *           = sum{j=0->inf} A(j) h^j sum{i=0->inf} B(i) h^i
 *
 *           = sum{k=0->inf} sum{j=0->k} A[j]B[k - j] h^k     where k = i + j  ==>  i = k - j
 *
 *  ==> C[k] = sum{j=0->k} A[j].B[k-j]     perhaps implemented by a static/private function cauchy(A, B, k, ...)
 *
 *  Dual numbers:
 *      C[0] = A[0].B[0]
 *      C[1] = A[0].B[1] + A[1].B[0]
 */

/*
 * Returns a pointer to kth element of the square of U, result stored and returned in variable S, NO JET STORAGE
 *
 *  S = U.U
 *
 *  S = sum{j=0->k} U[j].U[k-j]
 */
mpfr_t *t_sqr (series U, int k);

/*
 * Returns a pointer to kth element of the product of U and V, result stored in variable P, NO JET STORAGE
 *
 *  P = U.V
 *
 *  P = sum{j=0->k} U[j].V[k-j]
 */
mpfr_t *t_prod (series U, series V, int k);

/*
 * Returns a pointer to kth element of U / V, results accumulated in jet Q, DOMAIN RESTRICTION v[0] != 0.0
 *
 *     Q = U / V ==> U = Q.V
 *
 *                U[k] = sum{j=0->k} Q[j].V[k-j]
 *
 *                     = sum{j=0->k-1} Q[j].V[k-j] + Q[k].V[0]
 *
 *                Q[k] = (U[k] - sum{j=0->k-1} Q[j].V[k-j]) / V[0]
 *
 *                     =  U[0] / V[0]                                 if k == 0
 *
 *                     = (U[k] - sum{j=0->k-1} Q[j].V[k-j]) / V[0]    otherwise
 */
mpfr_t *t_quot (series Q, series U, series V, int k);

/*
 * Returns a pointer to kth element of 1 / V, results accumulated in jet I, DOMAIN RESTRICTION v[0] != 0.0
 *
 * from quotient, I[k] = 1.0 / V[0]                                   if k == 0
 *
 *                I[k] = - sum{j=0->k-1} I[j].V[k-j] / V[0]           otherwise
 */
mpfr_t *t_inv (series I, series V, int k);

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
 */
mpfr_t *t_sqrt (series r, series U, int k);

/*
 * Applying the chain rule for the derivative of a composed function f(u) creates another Cauchy product:
 *
 *           F' = (df/du).U' = H.U'    where H = df/du
 *
 *  Using F'[k] = (k+1) F[k+1]   (THE IDENTITY from earlier)
 *
 * ==>  F'[k-1] = k F[k], because we WANT F[k], and we can now replace F' with F, and U' with U as follows:
 *
 * Starting from the Cauchy product above, first rewrite it in terms of [k-1], then make the substitutions:
 *
 *        F'[k] = sum{j=0->k} H[j].U'[k-j]
 *
 *      F'[k-1] = sum{j=0->k-1} H[j].U'[k-1-j]
 *
 *        kF[k] = sum{j=0->k-1} H[j].(k-j)U[k-j]       if k > 0, need a mathematical function call for k == 0 (see Dual numbers)
 *
 *     ==> F[k] = sum{j=0->k-1} H[j].(k-j)U[k-j]/k     perhaps implemented by a static/private function d_cauchy(H, U, k, ...)
 *
 *  Dual numbers:
 *      F[0] = f(U[0])
 *      F[1] = H[0].U[1]
 */

/*
 * Returns a pointer to kth element of the exponential of U, results accumulated in jet E
 *
 *      E' = E.U'
 *
 *    E[k] = sum{j=0->k-1} E[j].(k-j)U[k-j]/k
 *
 *  Dual numbers:
 *      E[0] = exp(U[0])
 *      E[1] = E[0].U[1]
 */
mpfr_t *t_exp (series E, series U, int k);

/*
 * Returns a pair of pointers to kth elements of the sine and cosine of U, results accumulated in jets S and C
 *
 *      S' =       C.U'
 *      C' = (+/-) S.U'     + for cosh (g == HYP), - for cos (g == TRIG)
 *
 *    S[k] = sum{j=0->k-1}       C[j].(k-j)U[k-j]/k
 *    C[k] = sum{j=0->k-1} (+/-) S[j].(k-j)U[k-j]/k
 *
 *  Dual numbers:
 *      S[0] =      cos(U[0])
 *      C[0] = (+-) sin(U[0])
 *      S[1] =      C[0].U[1]
 *      C[1] = (+-) S[0].U[1]
 */
tuple t_sin_cos (series S, series C, series U, int k, geometry g);

/*
 * Returns a pair of pointers to kth elements of the tangent and squared secant of U, results accumulated in jets T and S2
 *
 *      T' =       S2.U'
 *     S2' = (+/-)2 T.T'    + for sec^2 (g == TRIG), - for sech^2 (g == HYP)
 *
 *    T[k] = sum{j=0->k-1}       S2[j].(k-j)U[k-j]/k
 *   S2[k] = sum{j=0->k-1} (+/-)2 T[j].(k-j)T[k-j]/k
 */
tuple t_tan_sec2 (series T, series S2, series U, int k, geometry g);

/*
 * Returns a pointer to kth element of P = U^a (where a is scalar), results accumulated in jet P, DOMAIN RESTRICTION U[0] > 0.0
 *
 *                                    P'= U^a' = a U^(a-1).U'
 *                                      U.U^a' = a U^a.U'
 *                                        U.P' = a P.U'
 *
 *              sum{j=0->k-1} U[j].(k-j)P[k-j] =  a sum{j=0->k-1} P[j].(k-j)U[k-j]
 *
 * U[0].kP[k] + sum{j=1->k-1} U[j].(k-j)P[k-j] =  a sum{j=0->k-1} P[j].(k-j)U[k-j]
 *
 *                                  U[0].kP[k] =  a sum{j=0->k-1} P[j].(k-j)U[k-j]   - sum{j=1->k-1} U[j].(k-j)P[k-j]
 *
 *                                        P[k] = (a sum{j=0->k-1} P[j].(k-j)U[k-j]/k - sum{j=1->k-1} U[j].(k-j)P[k-j]/k) / U[0]
 *
 *  Dual numbers:
 *      P[0] = pwr(U[0], a)
 *      P[1] = a * P[0].U[1] / U[0]
 */
mpfr_t *t_pwr (series P, series U, mpfr_t a, int k);

/*
 * Returns a pointer to kth element of the natural logarithm of U, result accumulated in jet L, DOMAIN RESTRICTION U[0] > 0.0
 *
 *                     L' = U' / U
 *                     U' = U.L'
 *
 *                 k.U[k] = sum{j=0->k-1} U[j].(k-j)L[k-j]
 *
 *                        = sum{j=1->k-1} U[j].(k-j)L[k-j] + U[0].kL[k]
 *
 *                   L[k] = (U[k] - sum{j=1->k-1} U[j].(k-j)L[k-j]/k) / U[0]
 */
mpfr_t *t_ln (series L, series U, int k);

