#!/bin/sh
#
#  Example: ./profile thomas 15 10 0.1 300000 1 0 0 .19
#
#  (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file

model=$1
shift
args=$*

progname="tsm_${model}_gprof"

gcc -Wall -pg -g -o $progname tsm-${model}.c taylor-ode.c main.c -lm

time -p ./$progname $args >/dev/null

gprof -p -b ./$progname gmon.out
