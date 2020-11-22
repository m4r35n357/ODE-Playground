/*
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

typedef long double real;

static int diverged (real x1, real y1, real z1, real x2, real y2, real z2, real threshold) {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)) > threshold;
}

int main(int argc, char **argv) {
    assert(argc >= 4);
    FILE *fileA, *fileB;
    if ((fileA = fopen(argv[1], "r")) == NULL || (fileB = fopen(argv[2], "r")) == NULL) {
        fprintf(stderr, "divergence: Cannot read data files!\n");
        exit(1);
    }
    real xA, yA, zA, xB, yB, zB, t;
    char bl[2];
    for (int i = 3; i < argc; i++) {
        real threshold = strtold(argv[i], NULL);
        while(fscanf(fileA, "%Le %Le %Le %Le %s %s %s\n", &xA, &yA, &zA, &t, bl, bl, bl) != EOF &&
              fscanf(fileB, "%Le %Le %Le %Le %s %s %s\n", &xB, &yB, &zB, &t, bl, bl, bl) != EOF) {
            if (diverged(xA, yA, zA, xB, yB, zB, threshold)) {
                fprintf(stdout, "%s %.1Le %s %.3Lf\n", "Threshold:", threshold, "t:", t);
                break;
            }
        }
    }
}
