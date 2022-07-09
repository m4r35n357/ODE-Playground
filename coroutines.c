/*
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int sq_gen (int a, int b);

int sq_gen (int low, int high) {
    static int i, resume = 0;
    if (resume) goto resume; else resume = 1;
    printf("Initializing, start = %d, end = %d\n", low, high);
    for (i = low; i <= high; i += 1) {
        printf("Looping, counter = %d\n", i);
        return i * i;
        resume:;
    }
    resume = 0;
    return 0;
}

int main(int argc, char **argv) {
    int square;

    assert(argc == 3);
    while ((square = sq_gen((int)strtol(argv[1], NULL, 10), (int)strtol(argv[2], NULL, 10)))) {
        printf("main() square = %d\n", square);
    }

    return 0;
}
