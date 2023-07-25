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
                "\x1B[1;31mFAIL\x1B[0;37m \x1B[1;37m%s\x1B[0;37m \x1B[0;35m%s\x1B[0;37m:\x1B[0;35m%i\x1B[0;37m\n", \
                #x, __FILE__, __LINE__); \
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

void tsm (int places, int n, mpfr_t h, int steps, mpfr_t x0, mpfr_t y0, mpfr_t z0, parameters *P, clock_t since);

/*
 * For returning x, y, z velocities from the model
 */
typedef struct {
    mpfr_t x;
    mpfr_t y;
    mpfr_t z;
} triplet;

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
 * Returns value if k is 0, and zero otherwise.  For handling _additive_ constants.
 */
mpfr_t *t_const (int n, mpfr_t a);

/*
 * Returns a pointer to kth element of the absolute value of U, no user-supplied jet storage needed
 */
mpfr_t *t_abs (series U, int k);

/*
 * Returns a pointer to kth element of the product of U and V, no user-supplied jet storage needed
 */
mpfr_t *t_mul (series U, series V, int k);

/*
 * Returns a pointer to kth element of the square of U, no user-supplied jet storage needed
 */
mpfr_t *t_sqr (series U, int k);

/*
 * Returns a pointer to kth element of U / V, results stored in user-supplied jet Q, DOMAIN RESTRICTION v[0] != 0.0
 */
mpfr_t *t_div (series Q, series U, series V, int k);

/*
 * Returns a pointer to kth element of the square root of U, results stored in user-supplied jet R, DOMAIN RESTRICTION U[0] > 0.0
 */
mpfr_t *t_sqrt (series R, series U, int k);

/*
 * Returns a pointer to kth element of the exponential of U, results stored in user-supplied jet E
 */
mpfr_t *t_exp (series E, series U, int k);

/*
 * For returning combined recurrence values
 */
typedef struct {
    mpfr_t *a;
    mpfr_t *b;
} pair;

/*
 * Returns struct of pointers to kth elements of both sine and cosine of U, results stored in user-supplied jets S and C
 */
pair t_sin_cos (series S, series C, series U, int k, bool trig);

/*
 * Returns struct of pointers to kth elements of both tangent and squared secant of U, results stored in user-supplied jets T and S2
 */
pair t_tan_sec2 (series T, series S2, series U, int k, bool trig);

/*
 * Returns a pointer to kth element of P = U^a (where a is scalar), results stored in user-supplied jet P, DOMAIN RESTRICTION U[0] > 0.0
 */
mpfr_t *t_pwr (series P, series U, mpfr_t a, int k);

/*
 * Returns a pointer to kth element of the natural logarithm of U, results stored in user-supplied jet L, DOMAIN RESTRICTION U[0] > 0.0
 */
mpfr_t *t_ln (series L, series U, int k);

/*
 * Returns kth elements of arcsin(h) of U and 1 / DF_DU, results stored in user-supplied jets As and DU_DF
 */
pair t_asin (series AS, series DU_DF, series U, int k, bool trig);

/*
 * Returns kth elements of arccos(h) of U and 1 / DF_DU, results stored in user-supplied jets As and DU_DF
 */
pair t_acos (series AC, series G, series U, int k, bool trig);

/*
 * Returns kth elements of arctan(h) of U and 1 / DF_DU, results stored in user-supplied jets As and DU_DF
 */
pair t_atan (series AT, series G, series U, int k, bool trig);
