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
 * Client model data
 */
typedef struct Parameters parameters;

/*
 * The numerical base for string IO conversions
 */
extern const int BASE;

/*
 * Global rounding strategy for MPFR
 */
extern const mpfr_rnd_t RND;

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

void tsm (int places, int n, mpfr_t h, int steps, series3 *jets, parameters *P, clock_t since);

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
 * Returns a pointer to kth element of the product of U and V, no user-supplied jet storage needed
 *
 *     P = U.V
 *
 *  P[k] = sum{j=0->k} U[j].V[k-j]
 */
mpfr_t *t_mul (series U, series V, int k);

/*
 * Returns a pointer to kth element of the square of U, no user-supplied jet storage needed
 *
 *     S = U.U
 *
 *  S[k] = sum{j=0->k} U[j].U[k-j]
 */
mpfr_t *t_sqr (series U, int k);

/*
 * Returns a pointer to kth element of U / V, results stored in user-supplied jet QUOT
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
 * If U == NULL, returns kth element of 1 / V, results stored in user-supplied jet Q,
 *
 * from above,    Q[k] = 1.0 / V[0]                                   if k == 0
 *
 *                Q[k] = - sum{j=0->k-1} Q[j].V[k-j] / V[0]           otherwise
 */
mpfr_t *t_div (series QUOT, series U, series V, int k);

/*
 * Returns a pointer to kth element of the square root of U, results stored in user-supplied jet ROOT
 *
 *     U = R.R
 *
 *  U[k] = sum{j=0->k} R[j].R[k - j]
 *
 *       = sum{j=1->k-1} R[j].R[k - j] + 2.R[k].R[0]
 *
 *  R[k] = (U[k] - sum{j=1->k-1} R[j].R[k-j]) / 2.R[0]
 */
mpfr_t *t_sqrt (series ROOT, series U, int k);

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
 *      SIN' =  COS.U' = COSH.U'
 *      COS' = -SIN.U' = SINH.U'
 *
 *    SIN[k] = sum{j=0->k-1}     COS[j].(k-j).U[k-j]/k
 *    COS[k] = sum{j=0->k-1} +/- SIN[j].(k-j).U[k-j]/k
 */
pair t_sin_cos (series SIN, series COS, series U, int k, bool trig);

/*
 * Returns pointers to kth elements of both tangent and squared secant of U, results stored in user-supplied jets TAN and SEC2
 *
 *      TAN' =  SEC^2.U' = SECH^2.U'
 *     SEC2' = 2.TAN.TAN' = -2.TANH.TANH'
 *
 *    TAN[k] = sum{j=0->k-1}     SEC2[j].(k-j).U[k-j]/k
 *   SEC2[k] = sum{j=0->k-1} +/-2 TAN[j].(k-j).TAN[k-j]/k
 */
pair t_tan_sec2 (series TAN, series SEC2, series U, int k, bool trig);

/*
 * Returns a pointer to kth element of the natural logarithm of U, results stored in user-supplied jet LN
 *
 * This is the simplest "reverse" example, where the wanted quantity is now on the RHS, and G = DU_DF!
 *
 *         LN' = [1 / U].U'
 *          U' = U.LN'
 *
 *       LN[k] = (U[k] - sum{j=1->k-1} U[j].(k-j).LN[k-j]/k) / U[0]     "reverse"
 */
mpfr_t *t_ln (series LN, series U, int k);

/*
 * Returns pointers to kth elements of arcsin/arsinh of SIN, and G == DF_DU, results stored in user-supplied jets U and G
 *
 *        SIN' =     G.U'  ==>  G = COS(U) = COSH(U)
 *          G' =  -SIN.U' = SINH.U'
 *
 *        U[k] = (SIN[k] - sum{j=1->k-1} G[j].(k-j).U[k-j]/k) / G[0]
 *
 *        G[k] = sum{j=0->k-1} +/- SIN[j].(k-j).U[k-j]/k
 */
pair t_asin (series U, series G, series SIN, int k, bool trig);

/*
 * Returns pointers to kth elements of arccos/arcosh of COS, and G == DF_DU, results stored in user-supplied jets U and G
 *
 *        COS' =   G.U'  ==>  G = -SIN(U) = SINH(U)
 *          G' = COS.U' = COSH.U'
 *
 *        U[k] = (COS[k] -/+ sum{j=1->k-1} G[j].(k-j).U[k-j]/k) / G[0]
 *
 *        G[k] = sum{j=0->k-1} COS[j].(k-j).U[k-j]/k
 */
pair t_acos (series U, series G, series COS, int k, bool trig);

/*
 * Returns pointers to kth elements of arctan/artanh of TAN, and G == DF_DU, results stored in user-supplied jets U and G
 *
 *        TAN' =     G.U'  ==>  G = SEC^2(U) = SECH^2(U)
 *          G' = 2.TAN.TAN' = -2.TANH.TANH'
 *
 *        U[k] = (TAN[k] - sum{j=1->k-1} G[j].(k-j).U[k-j]/k) / G[0]
 *
 *        G[k] = sum{j=0->k-1} +/-2 TAN[j].(k-j).TAN[k-j]/k
 */
pair t_atan (series U, series G, series TAN, int k, bool trig);

/*
 * Returns kth element of P = U^a (where a is scalar), results stored in user-supplied jet PWR
 *
 *        dp/dt =   (dp/du).(du/dt)
 *   PWR'= U^a' = a.U^(a-1).U'
 *       U.U^a' = a.U^a.U'
 *       U.PWR' = a.PWR.U'
 *
 *          P_k = a * sum{j=0->k-1} PWR[j].(k-j).U[k-j]/k                   (RHS)
 *
 *       PWR[k] = (P_k - sum{j=1->k-1} U[j].(k-j).PWR[k-j]/k) / U[0]        (LHS)
 *
 *              = (a * sum{j=0->k-1} PWR[j].(k-j).U[k-j] - sum{j=1->k-1} U[j].(k-j).PWR[k-j]) / k.U[0]
 */
mpfr_t *t_pwr (series PWR, series U, mpfr_t a, int k);
