#!/bin/sh

PROGNAME='tsm_lorenz_gprof'

gcc -Wall -pg -g tsm-lorenz.c taylor-ode.c main.c -o $PROGNAME -lm

./$PROGNAME 15 10 .01 1000000 -15.8 -17.48 35.64 10 28 8 3 >/dev/null

gprof -p -b ./$PROGNAME gmon.out
