/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include "real.h"

/*
 * Diffentiate real numbers from jets
 */
typedef long double *series;

/*
 * Prints an index column, and x, y, z columns, into a single line
 */
void t_output (long dp, real x, real y, real z, real t, char *x_label, char *y_label, char *z_label) ;

/*
 * Bulk set ODE parameters from the command line arguments (5 onwards)
 */
void t_params (char **argv, int count, ...);

/*
 * Creates a jet of the specified size
 */
series t_jet (long size);

/*
 * Safely and efficiently evaluate a polynomial of length n, with the coefficients in S, and the variable in h
 */
real t_horner (series S, long n, real h);

/*
 * The Taylor Series Method (TSM) in brief; to solve a system of Ordinary Differential Equations defined by:
 *
 *                         x'(t) = ode(x(t))
 *
 * We plug x(t) into the equations, then somehow use the value of x'(t) to estimate the next x value in the sequence.
 * 
 * Now, evolution of a time-varying quantity can be represented as a Taylor Series:
 *
 *                     x(t0 + h) = sum{k=0->inf} X[k].h^k     where X[0] = x(t0), X[k] = (d/dt)^k x(t0) / k! and h = t - t0  (A)
 *
 * This calculation is best performed using Horner's method. Similarly, the time derivative of x(t) can be represented as:
 *
 *                         x'(t) = sum{k=0->inf} X'[k].h^k                                                                   (B)
 *
 * So, X'[k] is the result of evaluating the ODE equation, but with all variables expressed as Taylor Series, X[k].
 *
 *                              X'[k] = ODE(X[k])
 *
 * Furthermore, by differentiating (A) and equating it to (B) we obtain an additional and very useful IDENTITY:
 * 
 *                     d/dt x(t) = sum{k=1->inf} k.X[k].h^(k-1)
 *
 *                               = sum{k=0->inf} (k+1).X[k+1].h^k                                                            (C)
 *
 * Comparing (B) and (C),  X'[k] = (k+1).X[k+1]             *** this IDENTITY is also used in deriving the recurrences below ***
 *
 *                   ==>  X[k+1] = X'[k] / (k + 1)
 *
 * 1. Starting with initial values X[0], evaluate the ODE equations (X'[k]) using X[k], and recurrence relations where needed,
 *    then generate the next Taylor Series coefficient X[k+1] using the IDENTITY
 *
 * 2. Apply Horner's method to calculate the new values x(t0 + h), which become X[0] for the next time step.
 */

/*
 * Returns the value of a if k is 0, and zero otherwise.  For handling _additive_ constants.
 */
real t_const (real a, int k);

/*
 * Returns kth element of the absolute value of U, no user-supplied jet storage needed
 */
real t_abs (series U, int k);

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
 *  ==> C[k] = sum{j=0->k} A[j].B[k-j]     perhaps implemented by a static/private function cauchy(A, B, k, ...)
 */

/*
 * Returns kth element of the product of U and V, no user-supplied jet storage needed
 *
 *  P = U.V
 *
 *  P = sum{j=0->k} U[j].V[k-j]
 */
real t_prod (series U, series V, int k);

/*
 * Returns kth element of the square of U, no user-supplied jet storage needed
 *
 *  S = U.U
 *
 *  S = sum{j=0->k} U[j].U[k-j]
 */
real t_sqr (series U, int k);

/*
 * Returns kth element of U / V, results accumulated in jet Q, DOMAIN RESTRICTION v[0] != 0.0
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
real t_quot (series Q, series U, series V, int k);

/*
 * Returns kth element of 1 / V, results accumulated in jet I, DOMAIN RESTRICTION v[0] != 0.0
 *
 * from quotient, I[k] = 1.0 / V[0]                                   if k == 0
 *
 *                I[k] = - sum{j=0->k-1} I[j].V[k-j] / V[0]           otherwise
 */
real t_inv (series I, series V, int k);

/*
 * Returns kth element of the square root of U, results stored in user-supplied jet R, DOMAIN RESTRICTION U[0] > 0.0
*
 *    U = R.R
 *
 * U[k] = sum{j=0->k} R[j].R[k - j]
 *
 *      = sum{j=1->k-1} R[j].R[k - j] + 2.R[k].R[0]
 *
 * R[k] = (U[k] - sum{j=1->k-1} R[j].R[k-j]) / (2.R[0])
 */
real t_sqrt (series R, series U, int k);

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
 *        kF[k] = sum{j=0->k-1} dFdU[j].(k-j).U[k-j]         only if k > 0.  Use a mathematical function call for dFdU[0]
 *
 *     ==> F[k] = sum{j=0->k-1} dFdU[j].(k-j).U[k-j] / k     can be implemented by a static function f_k(dFdU, U, k, ...)
 */

/*
 * Returns kth element of the exponential of U, results stored in user-supplied jet E
 *
 *      E' = E.U'
 *
 *    E[k] = sum{j=0->k-1} E[j].(k-j).U[k-j]/k
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
 * Returns kth elements of both sine and cosine of U, results stored in user-supplied jets S and C
 *
 *      S' =       C.U'
 *      C' = (+/-) S.U'     + for cosh (g == HYP), - for cos (g == TRIG)
 *
 *    S[k] = sum{j=0->k-1}       C[j].(k-j).U[k-j]/k
 *    C[k] = sum{j=0->k-1} (+/-) S[j].(k-j).U[k-j]/k
 */
pair t_sin_cos (series S, series C, series U, int k, geometry g);

/*
 * Returns kth elements of both tangent and squared secant of U, results stored in user-supplied jets T and S2
 *
 *      T' =       S2.U'
 *     S2' = (+/-)2.T.T'    + for sec^2 (g == TRIG), - for sech^2 (g == HYP)
 *
 *    T[k] = sum{j=0->k-1}       S2[j].(k-j).U[k-j]/k
 *   S2[k] = sum{j=0->k-1} (+/-)2 T[j].(k-j).T[k-j]/k
 */
pair t_tan_sec2 (series T, series S2, series U, int k, geometry g);

/*
 * Returns kth element of P = U^a (where a is scalar), results stored in user-supplied jet P, DOMAIN RESTRICTION U[0] > 0.0
 *
 *                                      P'= U^a' = a.U^(a-1).U'
 *                                        U.U^a' = a.U^a.U'
 *                                          U.P' = a.P.U'
 *
 *               sum{j=0->k-1} U[j].(k-j).P[k-j] =  a.sum{j=0->k-1} P[j].(k-j).U[k-j]
 *
 * U[0].k.P[k] + sum{j=1->k-1} U[j].(k-j).P[k-j] =  a.sum{j=0->k-1} P[j].(k-j).U[k-j]
 *
 *                                   U[0].k.P[k] =  a.sum{j=0->k-1} P[j].(k-j).U[k-j]   - sum{j=1->k-1} U[j].(k-j).P[k-j]
 *
 *                                          P[k] = (a.sum{j=0->k-1} P[j].(k-j).U[k-j]/k - sum{j=1->k-1} U[j].(k-j).P[k-j]/k) / U[0]
 */
real t_pwr (series P, series U, real a, int k);

/*
 * Returns kth element of the natural logarithm of U, results stored in user-supplied jet L, DOMAIN RESTRICTION U[0] > 0.0
 *
 *                     L' = (1/U).U'
 *                     U' = U.L'
 *
 *                 k.U[k] = sum{j=0->k-1} U[j].(k-j).L[k-j]
 *
 *                        = sum{j=1->k-1} U[j].(k-j).L[k-j] + U[0].kL[k]
 *
 *                   L[k] = (U[k] - sum{j=1->k-1} U[j].(k-j).L[k-j]/k) / U[0]
 */
real t_ln (series L, series U, int k);

/*
 * For returning x, y, z values
 */
typedef struct {
    real x;
    real y;
    real z;
} components;

/*
 * For performing a simulation
 */
void tsm (int argc, char **argv, long dp, long n, real h, long steps, real x0, real y0, real z0);

/*
 * Obligatory client method signatures
 */
void *get_p (int argc, char **argv, long order);

components ode (series x, series y, series z, void *params, int k);
