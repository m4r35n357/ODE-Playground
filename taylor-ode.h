/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <time.h>

/*
 * The numerical base for string IO conversions
 */
extern const int BASE;

/*
 * Global rounding strategy for MPFR
 */
extern const mpfr_rnd_t RND;

/*
 * Pre-allocates some local MPFR variables
 */
void t_init (int dp);

/*
 * Type for Taylor Series coordinate jets
 */
typedef mpfr_t *series;

/*
 * Prints a line of data to stdout
 */
void t_output (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t h, int step, double cpu);

/*
 * Retrieves ODE parameters from the tail of the command (arguments 9 onwards)
 */
void t_params (char **argv, int count, ...);

/*
 * Cumulative CPU time in seconds from a clock() value
 */
double cpu (clock_t since);

/*
 * Creates a zeroed Taylor Series jet with the specified number of elements
 */
series t_jet (int size);

/*
 * Safely and efficiently evaluates a polynomial of degree n, with the coefficients in S, and the variable in h
 */
mpfr_t *t_horner (series S, int n, mpfr_t h);

/*
 * The Taylor Series Method (TSM) in brief; to solve a system of Ordinary Differential Equations defined by:
 *
 *                         v(t) = ode(x(t))                  (where x represents x, y, z in turn)
 *
 * We plug x(t) into the equations, then somehow use the value(s) of v(t) to estimate the next x value(s) in the sequence.
 *
 * Now, evolution of a time-varying quantity can be represented as a Taylor Series:
 *
 *             x(t0 + h) = x(t) = sum{k=0->inf} X[k].h^k    where X[0] = x(t0), X[k] = (d/dt)^k x(t0) / k! and h = t - t0  (A)
 *
 * This calculation is best performed using Horner's method. Similarly, the velocity can be represented as:
 *
 *                         v(t) = sum{k=0->inf} V[k].h^k                                                                   (B)
 *
 * Where V[k] is the result of evaluating the ODE equation, but with all variables expressed as Taylor Series, X[k].
 *
 *                         V[k] = ODE(X[k])
 *
 * Furthermore, by explicitly differentiating (A) wrt t, we obtain an alternative description of the velocity:
 *
 *                    d/dt x(t) = sum{k=1->inf} k.X[k].h^(k-1)
 *
 *                              = sum{k=0->inf} (k+1).X[k+1].h^k                                                           (C)
 *
 * Comparing (B) and (C),  V[k] = (k+1).X[k+1]        *** this IDENTITY is also used in deriving the recurrences below ***
 *
 *                   ==> X[k+1] = V[k] / (k + 1)
 *
 * 1. Starting with initial values X[0], evaluate the ODE equations (V[k]) using X[k], and recurrence relations where needed,
 *
 *    then generate the next Taylor Series coefficient X[k+1], up to X[n], using the IDENTITY
 *
 * 2. Apply Horner's method to calculate the new values x(t0 + h), which become X[0] for the next time step.
 */
void tsm (int argc, char **argv, int n, mpfr_t h, int steps, mpfr_t x0, mpfr_t y0, mpfr_t z0);

/*
 * For returning x, y, z velocities from the model
 */
typedef struct {
    mpfr_t x;
    mpfr_t y;
    mpfr_t z;
} components;

/*
 * Obligatory client method signatures
 */

/*
 * Get a blob of parameter data from the model to be passed back in from ode()
 */
void *get_p (int argc, char **argv, int order);

/*
 * Calculate the kth components of the velocity jet V, using the coordinate jets and the parameter data,
 *
 * together with the functions below as necessary.
 */
void ode (components *V, series X, series Y, series Z, void *P, int k);

/*
 * Basic Taylor Series functions
 */

/*
 * Returns value if k is 0, and zero otherwise.  For handling _additive_ constants.
 */
mpfr_t *t_const (mpfr_t value, int k);

/*
 * Returns a pointer to kth element of the absolute value of U, no user-supplied jet storage needed
 */
mpfr_t *t_abs (series U, int k);

/*
 * Taylor Series recurrence relationships
 */

/*
 * Cauchy product for C = A.B
 *
 *   if c(t) = a(t).b(t)
 *
 * then c(t) = sum{k=0->inf} C(k).h^k
 *
 *           = sum{j=0->inf} A(j).h^j sum{i=0->inf} B(i).h^i
 *
 *           = sum{k=0->inf} sum{j=0->k} A[j].B[k - j].h^k     where k = i + j  ==>  i = k - j
 *
 *  ==> C[k] = sum{j=0->k} A[j].B[k-j]
 */

/*
 * Returns a pointer to kth element of the product of U and V, no user-supplied jet storage needed
 *
 *  P = U.V
 *
 *  P = sum{j=0->k} U[j].V[k-j]
 */
mpfr_t *t_mul (series U, series V, int k);

/*
 * Returns a pointer to kth element of the square of U, no user-supplied jet storage needed
 *
 *  S = U.U
 *
 *  S = sum{j=0->k} U[j].U[k-j]
 */
mpfr_t *t_sqr (series U, int k);

/*
 * Returns a pointer to kth element of U / V, results stored in user-supplied jet Q, DOMAIN RESTRICTION v[0] != 0.0
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
 *
 * If U == NULL, returns a pointer to kth element of 1 / V, results stored in user-supplied jet I, DOMAIN RESTRICTION v[0] != 0.0
 *
 * from quotient, I[k] = 1.0 / V[0]                                   if k == 0
 *
 *                I[k] = - sum{j=0->k-1} I[j].V[k-j] / V[0]           otherwise
 */
mpfr_t *t_div (series Q, series U, series V, int k);

/*
 * Returns a pointer to kth element of the square root of U, results stored in user-supplied jet R, DOMAIN RESTRICTION U[0] > 0.0
 *
 *    U = R.R
 *
 * U[k] = sum{j=0->k} R[j].R[k - j]
 *
 *      = sum{j=1->k-1} R[j].R[k - j] + 2.R[k].R[0]
 *
 * R[k] = (U[k] - sum{j=1->k-1} R[j].R[k-j]) / (2.R[0])
 */
mpfr_t *t_sqrt (series R, series U, int k);

/*
 * Applying the chain rule for the derivative of a composed function f(u) creates another Cauchy product:
 *
 *           F' = (df/du).U' = dFdU.U'
 *
 *  Using F'[k] = (k+1).F[k+1]   (the IDENTITY from earlier)
 *
 * ==>  F'[k-1] = k F[k], because we WANT F[k], and we can now replace F' with F, and U' with U as follows:
 *
 * Starting from the Cauchy product:
 *
 *        F'[k] = sum{j=0->k} dFdU[j].U'[k-j]
 *
 *  rewrite it in terms of [k-1]:
 *
 *      F'[k-1] = sum{j=0->k-1} dFdU[j].U'[k-1-j]            obviously does not work for k = 0 !
 *
 * then make the IDENTITY substitutions
 *
 *        kF[k] = sum{j=0->k-1} dFdU[j].(k-j).U[k-j]         only for k > 0.  Use a mathematical function call for dFdU[0]
 *
 *     ==> F[k] = sum{j=0->k-1} dFdU[j].(k-j).U[k-j] / k
 */

/*
 * Returns a pointer to kth element of the exponential of U, results stored in user-supplied jet E
 *
 *      E' = E.U'
 *
 *    E[k] = sum{j=0->k-1} E[j].(k-j).U[k-j]/k
 */
mpfr_t *t_exp (series E, series U, int k);

/*
 * Selects either a trigonometric or hyperbolic version of the function
 */
typedef enum {TRIG, HYP} geometry;

/*
 * For returning combined recurrence values
 */
typedef struct {
    mpfr_t *a;
    mpfr_t *b;
} pair;

/*
 * Returns struct of pointers to kth elements of both sine and cosine of U, results stored in user-supplied jets S and C
 *
 *      S' =       C.U'
 *      C' = (+/-) S.U'     + for cosh (g == HYP), - for cos (g == TRIG)
 *
 *    S[k] = sum{j=0->k-1}       C[j].(k-j).U[k-j]/k
 *    C[k] = sum{j=0->k-1} (+/-) S[j].(k-j).U[k-j]/k
 */
pair t_sin_cos (series S, series C, series U, int k, geometry g);

/*
 * Returns struct of pointers to kth elements of both tangent and squared secant of U, results stored in user-supplied jets T and S2
 *
 *      T' =       S2.U'
 *     S2' = (+/-)2.T.T'    + for sec^2 (g == TRIG), - for sech^2 (g == HYP)
 *
 *    T[k] = sum{j=0->k-1}       S2[j].(k-j).U[k-j]/k
 *   S2[k] = sum{j=0->k-1} (+/-)2 T[j].(k-j).T[k-j]/k
 */
pair t_tan_sec2 (series T, series S2, series U, int k, geometry g);

/*
 * Returns a pointer to kth element of P = U^a (where a is scalar), results stored in user-supplied jet P, DOMAIN RESTRICTION U[0] > 0.0
 *
 *                                      P'= U^a' = a.U^(a-1).U'
 *                                        U.U^a' = a.U^a.U'
 *                                          U.P' = a.P.U'
 *
 *               sum{j=0->k-1} U[j].(k-j).P[k-j] =  a.sum{j=0->k-1} P[j].(k-j).U[k-j]
 *
 * U[0].k.P[k] + sum{j=1->k-1} U[j].(k-j).P[k-j] =  a.sum{j=0->k-1} P[j].(k-j).U[k-j]
 *
 *                                          P[k] = (a.sum{j=0->k-1} P[j].(k-j).U[k-j] - sum{j=1->k-1} U[j].(k-j).P[k-j]) / k.U[0]
 *
 *                                               = (a.sum{j=0->k-1} P[j].(k-j).U[k-j] - sum{j=1->k-1} j.P[j].U[k-j]) / k.U[0]        (by symmetry)
 *
 *                                               = sum{j=0->k-1} (a.(k-j)-j).P[j].U[k-j] / k.U[0]    ("extra" term is 0.P[0].U[k])
 */
mpfr_t *t_pwr (series P, series U, mpfr_t a, int k);

mpfr_t *t_ipwr (series p, series u, int a, int k);

/*
 * Returns a pointer to kth element of the natural logarithm of U, results stored in user-supplied jet L, DOMAIN RESTRICTION U[0] > 0.0
 *
 *                     L' = (1/U).U'
 *                     U' = U.L'
 *
 *                 k.U[k] = sum{j=0->k-1} U[j].(k-j).L[k-j]
 *
 *                        = sum{j=1->k-1} U[j].(k-j).L[k-j] + U[0].k.L[k]
 *
 *                   L[k] = (U[k] - sum{j=1->k-1} U[j].(k-j).L[k-j]/k) / U[0]
 *
 *                        = (U[k] - sum{j=1->k-1} j.L[j].U[k-j]/k) / U[0]                    (by symmetry)
 */
mpfr_t *t_ln (series L, series U, int k);

/*
 * Returns kth elements of arcsin(h) of U and 1 / DF_DU, results stored in user-supplied jets As and DU_DF
 *
 *       df/du = 1 / sqrt(1 +- U^2)
 *       du/df = sqrt(1 +- U^2)
 *
 *         AS' =   dUdF.U'
 *      du/df' = (+/-)U.AS'    - for arcsin (g == TRIG), + for arcsinh (g == HYP)
 *
 *       AS[k] = (U[k] - sum{j=1->k-1} j.AS[j].dUdF[k-j]/k) / dUdF[0]                        (by symmetry)
 *
 *     dUdF[k] = sum{j=0->k-1} U[j].(k-j).AS[k-j]/k
 */
pair t_asin (series AS, series DU_DF, series U, int k, geometry g);

/*
 * Returns kth elements of arccos(h) of U and 1 / DF_DU, results stored in user-supplied jets As and DU_DF
 *
 *       df/du = -1 / sqrt(1 - U^2) for arccos (g == TRIG), 1 / sqrt(u^2 - 1) for arccosh (g == HYP)
 *       du/df =    - sqrt(1 - U^2) for arccos (g == TRIG),     sqrt(u^2 - 1) for arccosh (g == HYP)
 *
 *         AC' =   dUdF.U'
 *      du/df' = (+/-)U.AC'    - for arccos (g == TRIG), + for arccosh (g == HYP)
 *
 *       AC[k] = (U[k] - sum{j=1->k-1} j.AC[j].dUdF[k-j]/k) / dUdF[0]                        (by symmetry)
 *
 *     dUdF[k] = sum{j=0->k-1} U[j].(k-j).AC[k-j]/k
 */
pair t_acos (series AC, series G, series U, int k, geometry g);

/*
 * Returns kth elements of arctan(h) of U and 1 / DF_DU, results stored in user-supplied jets As and DU_DF
 *
 *       df/du = 1 / (1 +- U^2)
 *       du/df = (1 +- U^2)
 *
 *         AT' =     dUdF.U'
 *      du/df' = (+/-)2.U.U'    + for arctan (g == TRIG), - for arctanh (g == HYP)
 *
 *       AT[k] = (U[k] - sum{j=1->k-1} j.AT[j].dUdF[k-j]/k) / dUdF[0]                        (by symmetry)
 *
 *     dUdF[k] = sum{j=0->k-1} (+/-)2 U[j].(k-j).U[k-j]/k
 */
pair t_atan (series AT, series G, series U, int k, geometry g);
