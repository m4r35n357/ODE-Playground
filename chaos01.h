/*
 * Automatic search for chaos in ODE solutions - Interface
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

typedef enum {NO=0, YES=1} choice;

typedef enum {PC=0, MSD=1, KVC=2, K=3} function;

const long double PI = 3.14159265358979323846;

long double random_pi (void);

long double mean (long n, long double *data);

long double cov (long n, long double *a, long double *b);

long double corr (long n, long double *a, long double *b);

int compare (const void *a, const void *b);

long double median (int n, long double *array);

void import_data (long *n, long double *data, long column);

void translation_variables (long double c, long n, long double *data, long double *pc, long double *qc, choice print);

void mean_square_displacement (long double c, long n, long double *data, long double *pc, long double *qc,
                                long *n_cut, long double *mc, long double *dc, long double *xi, choice print);

