/*
 * Interface for solving systems of Ordinary Differential Equations
 * using the Taylor Series Method (TSM)
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <time.h>
#include <stdbool.h>

/*
 * Report program arguments in colour
 */
#define PRINT_ARGS(argc, argv) \
fprintf(stderr, "argc: \x1B[1;37m%d\x1B[0;37m, argv: [ \x1B[0;36m", argc); \
for (int i = 0; i < argc; i++) { \
    fprintf(stderr, "%s ", argv[i]); \
} \
fprintf(stderr, "\x1B[0;37m]\n");

/*
 * Unavoidable assert(), in colour
 */
#define CHECK(x) \
do { \
    if(!(x)) { \
        fprintf(stderr, \
            "\x1B[1;31mFAILED\x1B[0;37m \x1B[1;37m%s\x1B[0;37m %s() \x1B[1;37m%s\x1B[0;37m:\x1B[1;37m%i\x1B[0;37m\n", \
            #x, __func__, __FILE__, __LINE__); \
        exit(1); \
    } \
} while (0)

/*
 * Number base for integer conversions
 */
#define BASE 10

/*
 * Default MPFR rounding mode
 */
#define RND MPFR_RNDN

/*
 * Client model data
 */
typedef struct Parameters parameters;

/*
 * Type for Taylor Series coordinate jets
 */
typedef mpfr_t *series;

/*
 * Combined x, y, z series
 */
typedef struct triple_s {
    series x, y, z;
} series3;

/*
 * For returning x, y, z velocities from the model
 */
typedef struct triple_m {
    mpfr_t x, y, z;
} triplet;

/*
 * For returning combined recurrence values
 */
typedef struct pair_m {
    mpfr_t *a, *b;
} pair;

/*
 * The Taylor Series Method (TSM) in brief; to solve a system of Ordinary Differential Equations defined by:
 *
 *                         v(t) = ode(x(t))                  (where x represents x, y, z coordinates)
 *
 * We plug x(t) into the equations, then use the value(s) of the first derivatives v(t) to estimate the next x value(s) in the sequence.
 *
 * Now, a time-varying quantity x can be represented locally (instantaneously) as a Taylor Series:
 *
 *             x(t0 + h) = x(t) = sum{k=0->inf} X[k].h^k    where X[0] = x(t0), X[k] = (d/dt)^k x(t0) / k! and h = t - t0  (A)
 *
 * (this summation is best performed using Horner's method.) Similarly, the velocity v can be represented as:
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
 * which (together with all lower X derivatives) we then use to find V[k+1], and so on . . ..
 */

/*
 * Prints a line of data to stdout
 */
void t_out (mpfr_t x, mpfr_t y, mpfr_t z, mpfr_t h, int step, clock_t since);

/*
 * Retrieves ODE parameters from the tail of the command (arguments 9 onwards)
 */
void t_params (char **argv, int count, ...);

/*
 * Creates a zeroed Taylor Series jet with the specified number of elements
 */
series t_jet (int size);

/*
 * Safely and efficiently evaluates a polynomial of degree n, with the coefficients in S, and the variable in h
 */
mpfr_t *t_horner (series S, int n, mpfr_t h);

/*
 * Initialize constants
 */
void tsm_init(int display_precision);

/*
 * Execute the Taylor Series Method
 */
void tsm (int n, mpfr_t h, int steps, series3 *jets, parameters *P, clock_t since);

/*
 * Get a blob of parameter data from the model to be passed into ode()
 */
parameters *get_p (int argc, char **argv, int order);

/*
 * Calculate the kth components of the velocity jet V, using the coordinate jets and the parameter data,
 * together with the functions below as necessary.
 */
void ode (triplet *V, series X, series Y, series Z, parameters *P, int k);

/*
 * Basic Taylor Series functions
 */

/*
 * Returns value if k is 0, and zero otherwise.  For handling _additive_ constants in ODE models.
 */
mpfr_t *t_const (int n, mpfr_t a);

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
 *  ==> C[k] = sum{j=0->k} A[j].B[k-j] = sum{j=0->k} A[k-j].B[j]
 */

/*
 * Returns a pointer to kth element of the product of U and V, no user-supplied jet storage needed
 *
 *     PROD = U.V
 *
 *   PROD_k = sum{j=0->k} U[j].V[k-j]
 */
mpfr_t *t_mul (series U, series V, int k);

/*
 * Returns a pointer to kth element of the square of U, no user-supplied jet storage needed
 *
 *     SQR = U.U
 *
 *   SQR_k = sum{j=0->k} U[j].U[k-j]
 */
mpfr_t *t_sqr (series U, int k);

/*
 * Returns a pointer to kth element of U / V, results stored in user-supplied jet QUOT
 *
 *     QUOT = U / V ==> U = QUOT.V
 *
 *                   U[k] = sum{j=0->k} QUOT[j].V[k-j]
 *
 *                        = sum{j=0->k-1} QUOT[j].V[k-j] + QUOT[k].V[0]
 *
 *                QUOT[k] = (U[k] - sum{j=0->k-1} QUOT[j].V[k-j]) / V[0]
 *
 *                        =  U[0] / V[0]                                    if k == 0
 *
 *                        = (U[k] - sum{j=0->k-1} QUOT[j].V[k-j]) / V[0]    otherwise
 *
 * If U == NULL, returns kth element of 1 / V, results stored in user-supplied jet QUOT,
 *
 * from above,    QUOT[k] = 1.0 / V[0]                                      if k == 0
 *
 *                QUOT[k] = - sum{j=0->k-1} QUOT[j].V[k-j] / V[0]           otherwise
 */
mpfr_t *t_div (series QUOT, series U, series V, int k);

/*
 * Returns a pointer to kth element of the square root of U, results stored in user-supplied jet ROOT
 *
 *        U = ROOT.ROOT
 *
 *     U[k] = sum{j=0->k} ROOT[j].ROOT[k - j]
 *
 *          = sum{j=1->k-1} ROOT[j].ROOT[k - j] + 2.ROOT[k].ROOT[0]
 *
 *  ROOT[k] = (U[k] - sum{j=1->k-1} ROOT[j].ROOT[k-j]) / 2.ROOT[0]
 */
mpfr_t *t_sqrt (series ROOT, series U, int k);

/*
 * Applying the chain rule for the derivative of a composed function f(u(t)) creates another Cauchy product:
 *
 *           F' = (df/dt) = (df/du).(du/dt) = dFdU.U' = G.U'
 *
 *  Using F'[k] = (k+1).F[k+1]   (the IDENTITY from earlier)
 *
 * ==>  F'[k-1] = k.F[k], because we WANT F[k], and we can now replace F' with F, and U' with U as follows:
 *
 * Starting from the Cauchy product:
 *
 *        F'[k] = sum{j=1->k} G[j].U'[k-j]              differentiated series has one fewer terms
 *
 *  rewrite it with j index starting from 0, for "neatness":
 *
 *      F'[k-1] = sum{j=0->k-1} G[j].U'[k-1-j]
 *
 * then make the IDENTITY substitutions
 *
 *       k.F[k] = sum{j=0->k-1} G[j].(k-j).U[k-j]         only for k > 0.  Use a mathematical function call for G[0]
 *
 *     ==> F[k] = sum{j=0->k-1} G[j].(k-j).U[k-j] / k                             "forward"
 *
 *              = G[0].U[k] + sum{j=1->k-1} G[j].(k-j).U[k-j] / k
 *
 *     ==> U[k] = (F[k] - sum{j=1->k-1} G[j].(k-j).U[k-j]) / k) / dFdU[0]         "inverse"
 */

/*
 * Returns a pointer to kth element of the exponential of U, results stored in user-supplied jet EXP
 *
 *      EXP' = EXP.U'
 *
 *    EXP[k] = sum{j=0->k-1} EXP[j].(k-j).U[k-j]/k
 */
mpfr_t *t_exp (series EXP, series U, int k);

/*
 * Returns pointers to kth elements of both sine and cosine of U, results stored in user-supplied jets SIN and COS
 *
 *        SIN' =     G.U'  ==>  G = COS or COSH
 *          G' =  -SIN.U' or SINH.U'
 *
 *     SINH[0] = sinh(U[0])
 *
 *     COSH[0] = cosh(U[0])
 *
 *    SINH[k] = sum{j=0->k-1} COSH[j].(k-j).U[k-j] / k
 *
 *    COSH[k] = sum{j=0->k-1} SINH[j].(k-j).U[k-j] / k
 */
pair t_sin_cos (series SIN, series COS, series U, int k, bool trig);

/*
 * Returns pointers to kth elements of both tan(h) and sec(h)^2 of U, results stored in user-supplied jets TAN(H) and SEC(H)2
 *
 *      TAN' =  SEC^2.U' = SECH^2.U'
 *     SEC2' = 2.TAN.TAN' = -2.TANH.TANH'
 *
 *    TAN[0] = tan(U[0])
 *
 *   SEC2[0] = sec^2(U[0])
 *
 *    TAN[k] = sum{j=0->k-1}  SEC2[j].(k-j).U[k-j] / k
 *
 *   SEC2[k] = sum{j=0->k-1} 2.TAN[j].(k-j).TAN[k-j] / k
 */
pair t_tan_sec2 (series TAN, series SEC2, series U, int k, bool trig);

/*
 * Returns a pointer to kth element of the logarithm (inverse of EXP), results stored in user-supplied jet U
 *
 * This is the simplest "reverse" example, where the wanted quantity is now on the RHS
 *
 *        EXP' = EXP.U'  ==>  G = EXP
 *
 *        U[0] = ln(EXP[0])
 *
 *        U[k] = (EXP[k] - sum{j=1->k-1} EXP[j].(k-j).U[k-j] / k) / EXP[0]       "reverse"
 */
mpfr_t *t_ln (series U, series EXP, int k);

/*
 * Returns pointers to kth elements of arcsin/arsinh of SIN, and G == DF_DU, results stored in user-supplied jets U and G
 *
 *        SIN' =     G.U'  ==>  G = COS or COSH
 *          G' =  -SIN.U' or SINH.U'
 *
 *        U[0] = arsinh(SINH[0])
 *
 *     COSH[0] = cosh(U[0])
 *
 *        U[k] = (SINH[k] - sum{j=1->k-1} COSH[j].(k-j).U[k-j] / k) / COSH[0]    "reverse"
 *
 *     COSH[k] = sum{j=0->k-1} SINH[j].(k-j).U[k-j] / k                          "forward"
 */
pair t_asin (series U, series G, series SIN, int k, bool trig);

/*
 * Returns pointers to kth elements of arccos/arcosh (inverse of COSH), results stored in user-supplied jets U and G
 *
 *        COS' =   G.U'  ==>  G = -SIN or SINH
 *          G' = COS.U' or COSH.U'
 *
 *        U[0] = arcosh(COSH[0])
 *
 *     SINH[0] = sinh(U[0])
 *
 *        U[k] = (COSH[k] - sum{j=1->k-1} SINH[j].(k-j).U[k-j] / k) / SINH[0]    "reverse"
 *
 *     SINH[k] = sum{j=0->k-1} COSH[j].(k-j).U[k-j] / k                          "forward"
 */
pair t_acos (series U, series G, series COS, int k, bool trig);

/*
 * Returns pointers to kth elements of arctan/artanh (inverse of TAN), results stored in user-supplied jets U and G
 *
 *        TAN' =     G.U'  ==>  G = SEC^2 or SECH^2
 *          G' = 2.TAN.TAN' or -2.TANH.TANH'
 *
 *        U[0] = arctan(TAN[0])
 *
 *     SEC2[0] = sec^2(U[0])
 *
 *        U[k] = (TAN[k] - sum{j=1->k-1} SEC2[j].(k-j).U[k-j] / k) / SEC2[0]     "reverse"
 *
 *     SEC2[k] = sum{j=0->k-1} 2.TAN[j].(k-j).TAN[k-j] / k                       "forward"
 */
pair t_atan (series U, series G, series TAN, int k, bool trig);

/*
 * Returns kth element of P = U^a (where a is scalar), results stored in user-supplied jet PWR
 *
 *        dp/dt =   (dp/du).(du/dt)
 *   PWR'= U^a' = a.U^(a-1).U'
 *       U.U^a' = a.U^a.U'
 *       U.PWR' = a.PWR.U' = E_k
 *
 *          E_k = a * sum{j=0->k-1} PWR[j].(k-j).U[k-j]/k                   (RHS)
 *
 *       PWR[k] = (E_k - sum{j=1->k-1} U[j].(k-j).PWR[k-j]/k) / U[0]        (LHS)
 */
mpfr_t *t_pwr (series PWR, series U, mpfr_t a, int k);
