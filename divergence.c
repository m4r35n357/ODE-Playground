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
    real xA, yA, zA, xB, yB, zB, tA, tB, cpuA, cpuB;
    char *format = "%Le %Le %Le %Le %1s %1s %1s %Le", tag[2];
    for (int i = 3; i < argc; i++) {
        real threshold = strtold(argv[i], NULL);
        while(fscanf(fileA, format, &xA, &yA, &zA, &tA, tag, tag, tag, &cpuA) != EOF &&
              fscanf(fileB, format, &xB, &yB, &zB, &tB, tag, tag, tag, &cpuB) != EOF) {
            CHECK(tB == tA);
            if (sqrtl(SQR(xA - xB) + SQR(yA - yB) + SQR(zA - zB)) > threshold) {
                printf("%s %.1Le  %s %6.3Lf  %s %.3Lf\n", "threshold:", threshold, "t:", tB, "cpu:", cpuB);
                break;
            }
        }
    }
    return 0;
}
