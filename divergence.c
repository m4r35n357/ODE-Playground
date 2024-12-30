/*
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "real.h"

int main(int argc, char **argv) {
    PRINT_ARGS(argc, argv);
    CHECK(argc >= 4);
    FILE *fileA = fopen(argv[1], "r"); CHECK(fileA);
    FILE *fileB = fopen(argv[2], "r"); CHECK(fileB);
    long double xA, yA, zA, xB, yB, zB, tA, tB, cpuA, cpuB;
    char *format = "%Le %Le %Le %Le %Le";
    for (int i = 3; i < argc; i++) {
        long double threshold = strtold(argv[i], NULL);
        while(fscanf(fileA, format, &xA, &yA, &zA, &tA, &cpuA) != EOF &&
              fscanf(fileB, format, &xB, &yB, &zB, &tB, &cpuB) != EOF) {
            CHECK(tB == tA);
            if (sqrtl((xA - xB) * (xA - xB) + (yA - yB) * (yA - yB) + (zA - zB) * (zA - zB)) > threshold) {
                printf("%s %.1Le  %s %6.3Lf  %s %.3Lf\n", "threshold:", threshold, "t:", tB, "cpu:", cpuB);
                break;
            }
        }
    }
    return 0;
}
