/*
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

typedef long double real;

static real *get_array (long n) {
    return malloc(sizeof (real) * n);
}

static real separation (real x1, real y1, real z1, real x2, real y2, real z2) {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

static long read_data (FILE *fileA, FILE *fileB, real *xA, real *yA, real *zA, real *xB, real *yB, real *zB) {
    real tA, tB;
    long count = 0;
    while( fscanf(fileA, "%Le %Le %Le %Le\n", &xA[count], &yA[count], &zA[count], &tA) != EOF &&
           fscanf(fileB, "%Le %Le %Le %Le\n", &xB[count], &yB[count], &zB[count], &tB) != EOF) {
        count += 1;
    }
    return count;
}

int main(int argc, char **argv) {
    assert(argc = 3);
    long expected = strtol(argv[1], NULL, 10);
    real lower = strtold(argv[2], NULL);
    real upper = strtold(argv[3], NULL);
    real *xA = get_array(expected), *yA = get_array(expected), *zA = get_array(expected);
    real *xB = get_array(expected), *yB = get_array(expected), *zB = get_array(expected);
    FILE *fileA = fopen("/tmp/dataA", "r");
    FILE *fileB = fopen("/tmp/dataB", "r");
    if (fileA == NULL || fileB == NULL) { fprintf(stderr, "Cannot read files!\n"); exit(1); }
    long count = read_data(fileA, fileB, xA, yA, zA, xB, yB, zB);
    if (count == expected) {
        real max_s1 = -1.0;
        for (int i = 0; i < count; i++) {
            real s1 = separation(xA[i], yA[i], zA[i], xB[i], yB[i], zB[i]);
            max_s1 = s1 > max_s1 ? s1 : max_s1;
        }
        if (max_s1 <= lower) {
            fprintf(stdout, " LIMIT-CYCLE value = %.3Le\n", max_s1);
        } else if (max_s1 > upper) {
            fprintf(stdout, "     CHAOTIC value = %.3Le\n", max_s1);
        } else {
            fprintf(stdout, "UNCLASSIFIED value = %.3Le\n", max_s1);
        }
    } else if (count < expected) {
        fprintf(stdout, "   UNBOUNDED: value = -1.0\n");
    } else {
        fprintf(stdout, "INCORRECT DATA SIZE: read %ld lines, expected %ld\n", count, expected);
    }
}
