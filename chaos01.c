/*
 * Automatic search for chaos in ODE solutions
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "chaos01.h"

long double random_pi () {
    return ((long double)rand() / (long double)(RAND_MAX)) * PI;
}

long double mean (long n, long double data[]) {
    long double m = 0.0;
    for (int j = 1; j <= n; j++) {
        m += data[j];
    }
    return m / n;
}

long double cov (long n, long double a[], long double b[]) {
    long double a_bar = mean(n, a);
    long double b_bar = mean(n, b);
    long double c = 0.0;
    for (int j = 1; j <= n; j++) {
        c += (a[j] - a_bar) * (b[j] - b_bar);
    }
    return c / n;
}

long double corr (long n, long double *a, long double *b) {
    return cov(n, a, b) / sqrt(cov(n, a, a) * cov(n, b, b));
}

int compare (const void *a, const void *b) {
    long double fa = *(const long double *)a;
    long double fb = *(const long double *)b;
    return fa > fb ? 1 : -1;
}

long double median (int n, long double array[]) {
    qsort(array, n, sizeof(long double), compare);
    return n % 2 == 0 ? (array[n / 2 + 1] + array[n / 2]) / 2 : array[(n + 1) / 2];
}

void import_data (long *n, long double data[], long column) {
    long double x, y, z, t, *d;
    *n = 1;
    switch (column) {
        case 0 : d = &x; break;
        case 1 : d = &y; break;
        case 2 : d = &z; break;
        default :
            fprintf(stderr, ">>> ARG 2: invalid data column [%ld] - use one of [ 0 (x) | 1 (y) | 2 (z) ] <<<\n", column);
            exit(1);
    }
    while (scanf("%Le%Le%Le%Le", &x, &y, &z, &t) != EOF) {
        data[*n] = *d;
        *n += 1;
    }
    data[0] = NAN;  // temporary protection
}

void translation_variables (long double c, long n, long double data[], long double p[], long double q[], choice print) {
    for (int k = 1; k <= n; k++) {
        p[k] = 0.0;
        q[k] = 0.0;
        for (int j = 1; j <= k; j++) {
            p[k] += data[j] * cos(j * c);
            q[k] += data[j] * sin(j * c);
        }
        if (print) printf("%+.6d %+.12Le %+.12Le %+.12Le\n", k, data[k], p[k], q[k]);
    }
    p[0] = NAN; q[0] = NAN;  // temporary protection
}

void mean_square_displacement (long double c, long n, long *n_cut, long double data[], long double p[], long double q[],
                                long double m[], long double d[], long double xi[], choice print) {
    *n_cut = n / 10;
    long double e2 = mean(n, data) * mean(n, data);
    for (int k = 1; k <= *n_cut; k++) {
        m[k] = 0.0;
        for (int j = 1; j <= n - *n_cut; j++) {
            m[k] += ((p[j + k] - p[j]) * (p[j + k] - p[j]) + (q[j + k] - q[j]) * (q[j + k] - q[j])) / (n - *n_cut);
        }
        d[k] = m[k] - e2 * (1.0 - cos(k * c)) / (1.0 - cos(c));
        xi[k] = (long double)k;
        if (print) printf("%+.6d %+.12Le %+.12Le %+.12Le\n", k, m[k], d[k], xi[k]);
    }
    m[0] = NAN; d[0] = NAN; xi[0] = NAN;  // temporary protection
}

int main(int argc, char **argv) {
    long n, n_cut, nc = -1, command, column, random = NO;
    long double c = -1.0, data[100002], pc[100002], qc[100002], mc[10002], dc[10002], xi[10002],  kn[102], cn[102];

    assert(argc == 4 || argc == 5);
    command = strtol(argv[1], NULL, 10);
    column = strtol(argv[2], NULL, 10);
    switch (command) {
        case PC : case MSD :
            c = strtold(argv[3], NULL);
            if (c < 0.0 || c > PI) {
                fprintf(stderr, ">>> ARG 3: bad c value [%Lf] - should be 0.0 <= c <= PI <<<\n", c);
                exit(1);
            }
            break;
        case KVC : case K :
            nc = strtol(argv[3], NULL, 10);
            if (nc < 2 || nc % 2 == 1) {
                fprintf(stderr, ">>> ARG 3: bad c range [%ld] - should be even and greater than 2 <<<\n", nc);
                exit(1);
            }
            random = strtol(argv[4], NULL, 10);
            if (random != NO && random != YES) {
                fprintf(stderr, ">>> ARG 4: bad c domain [%ld] - use one of [ 0 (linear) | 1 (random) ] <<<\n", random);
                exit(1);
            }
            break;
        default :
            fprintf(stderr, ">>> ARG 1: command [%ld] not recognized - use one of [ 0 (PC) | 1 (MSD) | 2 (KVC) | 3 (K) ] <<<\n",
                    command);
            exit(1);
    }
    import_data(&n, data, column);
    switch (command) {
        case PC :
            translation_variables(c, n, data, pc, qc, YES);
            break;
        case MSD :
            translation_variables(c, n, data, pc, qc, NO);
            mean_square_displacement(c, n, &n_cut, data, pc, qc, mc, dc, xi, YES);
            break;
        case KVC :
            for (int p = 1; p <= nc; p++) {
                cn[p] = random == NO ? (p - 0.5) * PI / nc : random_pi();
                translation_variables(cn[p], n, data, pc, qc, NO);
                mean_square_displacement(cn[p], n, &n_cut, data, pc, qc, mc, dc, xi, NO);
                kn[p] = corr(n_cut, xi, dc);
                printf("%+.12Le %+.12Le\n", cn[p], kn[p]);
            }
            fprintf(stderr, "%+.12Le\n", median(nc, kn));
            break;
        case K :
            for (int r = 1; r <= nc; r++) {
                cn[r] = random == NO ? (r - 0.5) * PI / nc : random_pi();
                translation_variables(cn[r], n, data, pc, qc, NO);
                mean_square_displacement(cn[r], n, &n_cut, data, pc, qc, mc, dc, xi, NO);
                kn[r] = corr(n_cut, xi, dc);
                fprintf(stderr, "%+.12Le %+.12Le\n", cn[r], kn[r]);
            }
            printf("%+.12Le\n", median(nc, kn));
            break;
    }

    return 0;
}

