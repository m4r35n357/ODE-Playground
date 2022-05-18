/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 *
 * Example:  ./divergence /tmp/$USER/dataA /tmp/$USER/dataB 0.000000000001 0.000000001 0.000001 0.001 1.0
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "real.h"

int main(int argc, char **argv) {
    assert(argc >= 4);
    FILE *fileA, *fileB;
    if ((fileA = fopen(argv[1], "r")) == NULL || (fileB = fopen(argv[2], "r")) == NULL) {
        fprintf(stderr, "divergence: Cannot read data files!\n");
        exit(5);
    }
    real xA, yA, zA, xB, yB, zB, tA, tB, cpuA, cpuB;
    char bl[2];
    for (int i = 3; i < argc; i++) {
        real threshold = strtold(argv[i], NULL);
        while(fscanf(fileA, "%Le %Le %Le %Le %s %s %s %Le", &xA, &yA, &zA, &tA, bl, bl, bl, &cpuA) != EOF &&
              fscanf(fileB, "%Le %Le %Le %Le %s %s %s %Le", &xB, &yB, &zB, &tB, bl, bl, bl, &cpuB) != EOF) {
            assert(fabsl(tB - tA) <= 1.0e-99L);
            if (sqrtl((xA - xB) * (xA - xB) + (yA - yB) * (yA - yB) + (zA - zB) * (zA - zB)) > threshold) {
                fprintf(stdout, "%s %.1Le  %s %6.3Lf  %s %.3Lf\n", "threshold:", threshold, "t:", tB, "cpu:", cpuB);
                break;
            }
        }
    }
}
