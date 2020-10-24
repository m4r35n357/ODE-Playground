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

static long read_data (long expected, FILE *fileA, FILE *fileB, real *xA, real *yA, real *zA, real *xB, real *yB, real *zB) {
    real tA, tB;
    char bl[2];
    long count = 0;
    while(fscanf(fileA, "%Le %Le %Le %Le %s %s %s\n", &xA[count], &yA[count], &zA[count], &tA, bl, bl, bl) != EOF &&
          fscanf(fileB, "%Le %Le %Le %Le %s %s %s\n", &xB[count], &yB[count], &zB[count], &tB, bl, bl, bl) != EOF) {
        count += 1;
        if (count > expected) {
            fprintf(stderr, "chaos-distance-cns: Too much data!\n");
            exit(2);
        }
    }
    return count;
}

int main(int argc, char **argv) {
    assert(argc == 4);
    long expected = strtol(argv[1], NULL, 10);
    real *xA = get_array(expected), *yA = get_array(expected), *zA = get_array(expected);
    real *xB = get_array(expected), *yB = get_array(expected), *zB = get_array(expected);
    FILE *fileA = fopen("/tmp/dataA", "r");
    FILE *fileB = fopen("/tmp/dataB", "r");
    if (fileA == NULL || fileB == NULL) {
        fprintf(stderr, "chaos-distance-cns: Cannot read data files!\n");
        exit(1);
    }
    real lower = strtold(argv[2], NULL);
    real upper = strtold(argv[3], NULL);
    long count = read_data(expected, fileA, fileB, xA, yA, zA, xB, yB, zB);
    if (count == expected) {
        real s1 = -1.0, max_s1 = -1.0;
        for (int i = 0; i < count; i++) {
            s1 = separation(xA[i], yA[i], zA[i], xB[i], yB[i], zB[i]);
            max_s1 = s1 > max_s1 ? s1 : max_s1;
        }
        if (max_s1 <= lower) {
            fprintf(stdout, " %sLIMIT-CYCLE%s value = %.3Le classification = 0.0\n", "\x1b[1;32m", "\x1B[0m", max_s1);
        } else if (max_s1 > upper) {
            fprintf(stdout, "     %sCHAOTIC%s value = %.3Le classification = 1.0\n", "\x1b[1;31m", "\x1B[0m", max_s1);
        } else {
            fprintf(stdout, "%sUNCLASSIFIED%s value = %.3Le classification = 0.5\n", "\x1b[1;33m", "\x1B[0m", max_s1);
        }
    } else {
        fprintf(stdout, "   %sUNBOUNDED%s value = -1.0 classification = -0.1\n", "\x1b[1;37m", "\x1B[0m");
    }
}
